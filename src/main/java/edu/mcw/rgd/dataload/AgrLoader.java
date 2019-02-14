package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.FileDownloader;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by mtutaj on 3/14/2018.
 */
public class AgrLoader {

    XdbIdDAO xdao = new XdbIdDAO();
    OrthologRelationDao dao = new OrthologRelationDao();

    int invalidSpecies = 0;
    int noMatchHgnc = 0;
    int multiMatchHgnc = 0;
    int noMatchMgi = 0;
    int multiMatchMgi = 0;


    public static void main(String[] args) throws Exception {

        AgrLoader loader = new AgrLoader();
        loader.run();
    }

    void run() throws Exception {

        String[] taxons = { "Homo", "Mus", "Rattus", "Danio", "Drosophila", "Caenorhabditis", "Saccharomyces" };
        List<String> taxonList = new ArrayList<>(Arrays.asList(taxons));
        Collections.shuffle(taxonList);

        for( String taxon: taxonList ) {
            System.out.println("starting "+taxon);
            run(taxon);
        }
    }

    void run( String taxon ) throws Exception {

        final int rowsInBatch = 1000;
        String url = "https://www.alliancegenome.org/api/homologs/{{TAXON}}?stringencyFilter=stringent&rows="+rowsInBatch+"&start=";
        url = url.replace("{{TAXON}}", taxon);

        FileDownloader fd = new FileDownloader();
        fd.setDoNotUseHttpClient(true);

        for( int start=1; ; start += rowsInBatch ) {
            fd.setExternalFile(url+start);
            fd.setLocalFile("data/"+taxon+"_"+start+".json");
            String localFile = fd.downloadNew();
            int recordsParsed = handleFile(localFile);
            if( recordsParsed==0 ) {
                break;
            }
        }
    }

    int handleFile(String fileName) throws Exception {

        int linesInserted = 0;
        int linesUpdated = 0;

        JSONParser parser = new JSONParser();
        Object obj = parser.parse(new FileReader(fileName));
        JSONObject root = (JSONObject) obj;

        JSONArray records = (JSONArray) root.get("results");
        System.out.println("parsing "+fileName+"; records="+records.size());
        for( Object rec: records ) {
            JSONObject o = (JSONObject) rec;

            JSONObject gene1 = (JSONObject) o.get("gene");
            String speciesName = (String) gene1.get("speciesName");
            String geneSymbol = (String) gene1.get("symbol");
            String geneId = (String) gene1.get("geneID");
            Gene g1 = resolveGene(speciesName, geneSymbol, geneId);

            JSONObject gene2 = (JSONObject) o.get("homologGene");
            speciesName = (String) gene2.get("speciesName");
            geneSymbol = (String) gene2.get("symbol");
            geneId = (String) gene2.get("geneID");
            Gene g2 = resolveGene(speciesName, geneSymbol, geneId);

            boolean isBestRevScore = (boolean) o.get("best");
            boolean isBestScore = (boolean) o.get("bestReverse");

            String methodsMatched = "";
            JSONArray arr = (JSONArray) o.get("predictionMethodsMatched");
            for( Object item: arr ) {
                if( !methodsMatched.isEmpty() ) {
                    methodsMatched += ", ";
                }
                methodsMatched += item.toString();
            }

            String methodsNotCalled = "";
            arr = (JSONArray) o.get("predictionMethodsNotCalled");
            for( Object item: arr ) {
                if( !methodsNotCalled.isEmpty() ) {
                    methodsNotCalled += ", ";
                }
                methodsNotCalled += item.toString();
            }

            String methodsNotMatched = "";
            arr = (JSONArray) o.get("predictionMethodsNotMatched");
            for( Object item: arr ) {
                if( !methodsNotMatched.isEmpty() ) {
                    methodsNotMatched += ", ";
                }
                methodsNotMatched += item.toString();
            }

            // String confidence = (String) o.get("confidence");  ???
            String confidence = "stringent";

            int r = updateDb(g1.getRgdId(), g2.getRgdId(), confidence, isBestScore, isBestRevScore, methodsMatched, methodsNotMatched, methodsNotCalled);
            if( r==1 ) {
                linesInserted++;
            } else {
                linesUpdated++;
            }
        }
        return records.size();
    }

    Gene resolveGene(String speciesName, String geneSymbol, String geneId) throws Exception {

        int speciesTypeKey = SpeciesType.parse(speciesName);
        if( speciesTypeKey <= 0 ) {
            System.out.println("SP problem");
        }
        Gene gene = null;

        // first resolve by xdb id
        if( speciesTypeKey==SpeciesType.RAT ) {
            int geneRgdId = Integer.parseInt(geneId.substring(4));
            gene = dao.getGeneByRgdId(geneRgdId);
        } else {
            int xdbKey = speciesTypeKey==SpeciesType.HUMAN ? XdbId.XDB_KEY_HGNC
                    : speciesTypeKey==SpeciesType.MOUSE ? XdbId.XDB_KEY_MGD
                    : 63; // AGR_GENE

            List<Gene> genes = dao.getGenesByXdbId(geneId, xdbKey);
            if( genes.size()>1 ) {
                throw new Exception("multiple HGNC ids for "+geneId);
            }
            if( !genes.isEmpty() ) {
                gene = genes.get(0);
            }
        }

        // second, resolve by gene symbol
        if( gene==null ) {
            gene = dao.getGeneBySymbol(geneSymbol, speciesTypeKey);
        }

        // if gene still not found, insert it
        if( gene==null ) {
            gene = dao.insertAgrGene(speciesTypeKey, geneSymbol, geneId);
        } else {
            // fixup for zebrafish
            if( speciesTypeKey==SpeciesType.ZEBRAFISH ) {
                List<XdbId> xdbIds = xdao.getXdbIdsByRgdId(63, gene.getRgdId());
                if( xdbIds.isEmpty() ) {
                    XdbId xdbId = new XdbId();
                    xdbId.setAccId(geneId);
                    xdbId.setRgdId(gene.getRgdId());
                    xdbId.setXdbKey(63);
                    xdbId.setSrcPipeline("AgrOrtholog");
                    xdao.insertXdb(xdbId);
                }
            }
        }
        return gene;
    }

    void runOld() throws Exception {

        int linesInserted = 0;
        int linesUpdated = 0;

        String fname = "/rgd/orthology_RGD_1.0.0.0_1.json";
        //String fname = "/rgd/orthology_MGI_1.0.0.0_1.json";
        //String fname = "/rgd/orthology_Human_1.0.0.0_1.json";

        System.out.println("Loading file "+fname);

        JSONParser parser = new JSONParser();
        Object obj = parser.parse(new FileReader(fname));
        JSONObject root = (JSONObject) obj;

        JSONArray records = (JSONArray) root.get("data");
        System.out.println("parsing data; records="+records.size());
        for( Object rec: records ) {
            JSONObject o = (JSONObject) rec;
            String gene1 = (String) o.get("gene1");
            //String gene1Species = (String) o.get("gene1Species");
            String gene2 = (String) o.get("gene2");
            //String gene2Species = (String) o.get("gene2Species");

            // first species is always rat rgd id, f.e. "DRSC:RGD:1305631"
            int rgdId1 = getRgdId(gene1);
            if( rgdId1<=0 ) {
                continue;
            }

            // second species: we parse only human and mouse, f.e. "DRSC:HGNC:16486" "DRSC:MGI:1915283"
            int rgdId2 = getRgdId(gene2);
            if( rgdId2<=0 ) {
                continue;
            }

            String confidence = (String) o.get("confidence");
            boolean isBestRevScore = (boolean) o.get("isBestRevScore");
            boolean isBestScore = (boolean) o.get("isBestScore");

            String methodsMatched = "";
            JSONArray arr = (JSONArray) o.get("predictionMethodsMatched");
            for( Object item: arr ) {
                if( !methodsMatched.isEmpty() ) {
                    methodsMatched += ", ";
                }
                methodsMatched += item.toString();
            }

            String methodsNotCalled = "";
            arr = (JSONArray) o.get("predictionMethodsNotCalled");
            for( Object item: arr ) {
                if( !methodsNotCalled.isEmpty() ) {
                    methodsNotCalled += ", ";
                }
                methodsNotCalled += item.toString();
            }

            String methodsNotMatched = "";
            arr = (JSONArray) o.get("predictionMethodsNotMatched");
            for( Object item: arr ) {
                if( !methodsNotMatched.isEmpty() ) {
                    methodsNotMatched += ", ";
                }
                methodsNotMatched += item.toString();
            }

            int r = updateDb(rgdId1, rgdId2, confidence, isBestScore, isBestRevScore, methodsMatched, methodsNotMatched, methodsNotCalled);
            if( r==1 ) {
                linesInserted++;
            } else {
                linesUpdated++;
            }
        }

        System.out.println("===");
        System.out.println("lines inserted: "+linesInserted);
        System.out.println("lines updated: "+linesUpdated);
        System.out.println("lines with other species: "+invalidSpecies);
        System.out.println("lines HGNC no match: "+noMatchHgnc);
        System.out.println("lines HGNC multi match: "+multiMatchHgnc);
        System.out.println("lines MGI no match: "+noMatchMgi);
        System.out.println("lines MGI multi match: "+multiMatchMgi);

        System.out.println("===");
        System.out.println("run this: delete from agr_orthologs where last_update_date<sysdate-1");
    }

    int getRgdId(String geneId) throws Exception {
        // second species: we parse only human and mouse, f.e. "DRSC:HGNC:16486" "DRSC:MGI:1915283"
        int rgdId = 0;
        String accId;
        if( geneId.startsWith("DRSC:HGNC:") ) {
            accId = geneId.substring(5);

            List<Gene> genes = xdao.getActiveGenesByXdbId(XdbId.XDB_KEY_HGNC, accId);
            if( genes.isEmpty() ) {
                noMatchHgnc++;
                return -1;
            }
            if( genes.size()>1 ) {
                System.out.println("MULTI: "+accId);
                multiMatchHgnc++;
                return -2;
            }
            rgdId = genes.get(0).getRgdId();

        } else if( geneId.startsWith("DRSC:MGI:")) {
            accId = geneId.substring(5);

            List<Gene> genes = xdao.getActiveGenesByXdbId(XdbId.XDB_KEY_MGD, accId);
            if( genes.isEmpty() ) {
                noMatchMgi++;
                return -3;
            }
            if( genes.size()>1 ) {
                System.out.println("MULTI: "+accId);
                multiMatchMgi++;
                return -4;
            }
            rgdId = genes.get(0).getRgdId();

        } else if( geneId.startsWith("DRSC:RGD:")) {
            accId = geneId.substring(9);
            rgdId = Integer.parseInt(accId);
        } else {
            invalidSpecies++;
            return -5;
        }
        return rgdId;
    }

    int updateDb(int rgdId1, int rgdId2, String confidence, boolean isBestScore, boolean isBestRevScore,
                  String methodsMatched, String methodsNotMatched, String methodsNotCalled) throws Exception {

        String sql1 = "SELECT COUNT(0) FROM agr_orthologs WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";
        int dataExists =  xdao.getCount(sql1, rgdId1, rgdId2, methodsMatched);
        if( dataExists==0 ) {
            // insert
            final String sql = "INSERT INTO agr_orthologs (gene_rgd_id_1, gene_rgd_id_2, confidence, is_best_score, "+
                    "is_best_rev_score, methods_matched, methods_not_matched, methods_not_called, last_update_date) "+
                    "VALUES(?,?,?,?,?,?,?,?,SYSDATE)";
            xdao.update(sql, rgdId1, rgdId2, confidence, isBestScore?"Y":"N", isBestRevScore?"Y":"N",
                    methodsMatched, methodsNotMatched, methodsNotCalled);
            return 1;
        } else{
            final String sql = "UPDATE agr_orthologs SET confidence=?, is_best_score=?, is_best_rev_score=?, " +
                    "methods_not_matched=?, methods_not_called=?, last_update_date=SYSDATE " +
                    "WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";
            xdao.update(sql, confidence, isBestScore?"Y":"N", isBestRevScore?"Y":"N",
                    methodsNotMatched, methodsNotCalled, rgdId1, rgdId2, methodsMatched);
            return 0;
        }
    }
}
