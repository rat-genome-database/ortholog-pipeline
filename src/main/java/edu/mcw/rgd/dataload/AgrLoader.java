package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.object.MappingSqlQuery;

import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.text.ParseException;
import java.util.*;

/**
 * Created by mtutaj on 3/14/2018.
 */
public class AgrLoader {

    XdbIdDAO xdao = new XdbIdDAO();
    OrthologRelationDao dao = new OrthologRelationDao();

    Logger logDel = Logger.getLogger("deletedAgrOrthologs");
    Logger logIns = Logger.getLogger("insertedAgrOrthologs");

    int inserted = 0;
    private Set<String> processedSpecies;

    public void run() throws Exception {

        List<String> speciesList = new ArrayList<>(getProcessedSpecies());
        Collections.shuffle(speciesList);

        for( String speciesName: speciesList ) {
            System.out.println("starting "+speciesName);
            run(speciesName);
        }
    }

    void run( String speciesName ) throws Exception {

        // note: there could be disparity between oracle server time and the time on machine the pipeline is running,
        //  so to prevent deletion of valid orthologs (classified by the code as 'stale orthologs'),
        //  we use the as the cutoff timestamp the machine timestamp minus one day
        Date time0 = Utils.addDaysToDate(new Date(), -1); // time0 = current day-and-time less 1 day

        int speciesTypeKey = SpeciesType.parse(speciesName);
        String genus = SpeciesType.getOrganismGenus(speciesTypeKey); // f.e. 'Rattus'
        if( Utils.isStringEmpty(genus) ) {
            throw new Exception("ERROR: unknown species "+speciesName);
        }

        // delete json files from AGR that are older than 10 days
        int filesDeleted = deleteOldJsonFiles(genus, 10);
        if( filesDeleted>0 ) {
            System.out.println("deleted old json files: "+filesDeleted);
        }

        inserted = 0;

        final int rowsInBatch = 1000;
        String url = "https://www.alliancegenome.org/api/homologs/{{TAXON}}?stringencyFilter=stringent&rows="+rowsInBatch+"&start=";
        url = url.replace("{{TAXON}}", genus);

        FileDownloader fd = new FileDownloader();
        fd.setDoNotUseHttpClient(true);

        for( int start=1; ; start += rowsInBatch ) {
            fd.setExternalFile(url+start);
            fd.setLocalFile("data/"+genus+"_"+start+".json");
            String localFile = fd.downloadNew();
            int recordsParsed = handleFile(localFile, speciesTypeKey);
            if( recordsParsed==0 ) {
                break;
            }
        }

        System.out.println("inserted ortholog count for "+genus+": "+inserted);

        // delete stale orthologs
        String sql = "SELECT COUNT(0) FROM agr_orthologs WHERE "+
                " EXISTS( SELECT 1 FROM rgd_ids WHERE gene_rgd_id_1=rgd_id AND species_type_key=? )";
        int orthologCount = xdao.getCount(sql, speciesTypeKey);
        System.out.println("total ortholog count for "+genus+": "+orthologCount);


        sql = "SELECT * FROM agr_orthologs WHERE last_update_date<? "+
                "AND EXISTS( SELECT 1 FROM rgd_ids WHERE gene_rgd_id_1=rgd_id AND species_type_key=? )";
        MappingSqlQuery q = new MappingSqlQuery(xdao.getDataSource(), sql) {
            @Override
            protected Object mapRow(ResultSet rs, int rowNum) throws SQLException {
                int geneRgdId1 = rs.getInt("gene_rgd_id_1");
                int geneRgdId2 = rs.getInt("gene_rgd_id_2");
                String conf = rs.getString("confidence");
                String isBest = rs.getString("is_best_score");
                String isBestRev = rs.getString("is_best_rev_score");
                String methodsMatched = rs.getString("methods_matched");
                String methodsNotMatched = rs.getString("methods_not_matched");
                String methodsNotCalled = rs.getString("methods_not_called");
                Date createdDate = rs.getTimestamp("created_date");
                Date lastUpdateDate = rs.getTimestamp("last_update_date");

                logDel.debug(speciesName+"|RGD1:"+geneRgdId1+"|RGD2:"+geneRgdId2+"|"+conf+
                    "|IS_BEST:"+isBest+"|IS_BEST_REV:"+isBestRev+"|METHODS_MATCHED:"+methodsMatched+
                    "|METHODS_NOT_MATCHED:"+methodsNotMatched+"|METHODS_NOT_CALLED:"+methodsNotCalled+
                    "|CREATED:"+createdDate+"|LAST_UPDATED:"+lastUpdateDate);
                return null;
            }
        };
        q.declareParameter(new SqlParameter(Types.TIMESTAMP));
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.compile();
        List staleOrthologs = q.execute(time0, speciesTypeKey);
        System.out.println(" stale ortholog count for "+genus+": "+staleOrthologs.size());

        if( staleOrthologs.size() > orthologCount ) {
            System.out.println("*** WARN *** Cannot delete more than 10% of stale orthologs!");
        } else {
            sql = "DELETE FROM agr_orthologs WHERE last_update_date<? " +
                    "AND EXISTS( SELECT 1 FROM rgd_ids WHERE gene_rgd_id_1=rgd_id AND species_type_key=? )";
            int staleRowsDeleted = xdao.update(sql, time0, speciesTypeKey);
            System.out.println("stale rows deleted from AGR_ORTHOLOGS for " + genus + ": " + staleRowsDeleted);
        }
    }

    int handleFile(String fileName, int speciesTypeKey) throws Exception {

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
            String speciesName = (String) gene1.get("taxonId");
            String geneSymbol = (String) gene1.get("symbol");
            String geneId = (String) gene1.get("id");
            Gene g1 = resolveGene(speciesName, geneSymbol, geneId);
            if( g1==null ) {
                System.out.println("WARN: cannot resolve gene ["+geneSymbol+"] ["+geneId+"] "+speciesName);
                continue;
            }
            if( g1.getSpeciesTypeKey()!=speciesTypeKey ) {
                System.out.println("unexpected first species type key "+g1.getSpeciesTypeKey());
            }
            JSONObject gene2 = (JSONObject) o.get("homologGene");
            speciesName = (String) gene2.get("taxonId");
            geneSymbol = (String) gene2.get("symbol");
            geneId = (String) gene2.get("id");
            Gene g2 = resolveGene(speciesName, geneSymbol, geneId);
            if( g2==null ) {
                System.out.println("WARN: cannot resolve gene ["+geneSymbol+"] ["+geneId+"] "+speciesName);
                continue;
            }

            boolean isBestScore = o.get("best").equals("yes");
            boolean isBestRevScore = o.get("bestReverse").equals("Yes");

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
                inserted++;
            } else {
                linesUpdated++;
            }
        }
        return records.size();
    }

    Gene resolveGene(String speciesName, String geneSymbol, String geneId) throws Exception {

        // taxon fix up for species
        if( speciesName.equals("NCBITaxon:559292") ) {
            speciesName = "NCBITaxon:4932";
        }

        if( speciesName.startsWith("NCBITaxon:") ) {
            speciesName = "taxon:"+speciesName.substring(10);
        }

        int speciesTypeKey = SpeciesType.parse(speciesName);
        if( speciesTypeKey <= 0 ) {
            System.out.println("SP problem: speciesName ["+speciesName+" geneSymbol ["+geneSymbol+"] geneId ["+geneId+"]");
            return null;
        }
        Gene gene = null;

        // first resolve by xdb id
        int xdbKey = 63; // AGR_GENE
        List<Gene> genes = dao.getGenesByXdbId(geneId, xdbKey);
        if( genes.size()>1 ) {
            throw new Exception("multiple genes for "+geneId);
        }
        if( !genes.isEmpty() ) {
            gene = genes.get(0);
        }

        // second: special resolution for rat,mouse,human
        if( gene==null ) {
            if( speciesTypeKey==SpeciesType.RAT ) {
                int geneRgdId = Integer.parseInt(geneId.substring(4));
                gene = dao.getGeneByRgdId(geneRgdId);
            }
            else if( speciesTypeKey==SpeciesType.MOUSE ) {
                genes = dao.getGenesByXdbId(geneId, XdbId.XDB_KEY_MGD);
                if( !genes.isEmpty() ) {
                    gene = genes.get(0);
                }
            }
            else if( speciesTypeKey==SpeciesType.HUMAN ) {
                genes = dao.getGenesByXdbId(geneId, XdbId.XDB_KEY_HGNC);
                if( !genes.isEmpty() ) {
                    gene = genes.get(0);
                }
            }

            if( gene!=null ) {
                // gene resolved: insert AGR_GENE xdbid
                dao.insertAgrGeneXdbId(gene.getRgdId(), geneId);
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
                    dao.insertAgrGeneXdbId(gene.getRgdId(), geneId);
                }
            }
        }
        return gene;
    }

    int updateDb(int rgdId1, int rgdId2, String confidence, boolean isBestScore, boolean isBestRevScore,
                  String methodsMatched, String methodsNotMatched, String methodsNotCalled) throws Exception {

        String bestScore = isBestScore?"Y":"N";
        String bestRevScore = isBestRevScore?"Y":"N";

        String sql1 = "SELECT COUNT(0) FROM agr_orthologs WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";
        int dataExists =  xdao.getCount(sql1, rgdId1, rgdId2, methodsMatched);
        if( dataExists==0 ) {
            // insert
            final String sql = "INSERT INTO agr_orthologs (gene_rgd_id_1, gene_rgd_id_2, confidence, is_best_score, "+
                    "is_best_rev_score, methods_matched, methods_not_matched, methods_not_called, last_update_date) "+
                    "VALUES(?,?,?,?,?,?,?,?,SYSDATE)";
            xdao.update(sql, rgdId1, rgdId2, confidence, bestScore, bestRevScore,
                    methodsMatched, methodsNotMatched, methodsNotCalled);

            logIns.debug("RGD1:"+rgdId1+"|RGD2:"+rgdId2+"|"+confidence+
                    "|IS_BEST:"+bestScore+"|IS_BEST_REV:"+bestRevScore+"|METHODS_MATCHED:"+methodsMatched+
                    "|METHODS_NOT_MATCHED:"+methodsNotMatched+"|METHODS_NOT_CALLED:"+methodsNotCalled);
            return 1;
        } else{
            final String sql = "UPDATE agr_orthologs SET confidence=?, is_best_score=?, is_best_rev_score=?, " +
                    "methods_not_matched=?, methods_not_called=?, last_update_date=SYSDATE " +
                    "WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";

            xdao.update(sql, confidence, bestScore, bestRevScore,
                    methodsNotMatched, methodsNotCalled, rgdId1, rgdId2, methodsMatched);
            return 0;
        }
    }

    int deleteOldJsonFiles( String genus, int dateOffset ) throws ParseException {

        Date cutoffDate = Utils.addDaysToDate(new Date(), -dateOffset);
        int filesDeleted = 0;

        File dir = new File("data");
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith(genus) && name.endsWith(".json");
            }
        });

        for( File f: files ) {
            // delete files older than 10 days
            if( f.lastModified() < cutoffDate.getTime() ) {
                f.delete();
                filesDeleted++;
                continue;
            }
            // also delete files with 0 size
            if( f.length()==0 ) {
                f.delete();
                filesDeleted++;
            }
        }
        return filesDeleted;
    }

    public void setProcessedSpecies(Set<String> processedSpecies) {
        this.processedSpecies = processedSpecies;
    }

    public Set<String> getProcessedSpecies() {
        return processedSpecies;
    }
}
