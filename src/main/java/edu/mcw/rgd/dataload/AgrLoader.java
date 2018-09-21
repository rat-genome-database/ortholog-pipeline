package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.XdbId;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import java.io.FileReader;
import java.util.List;

/**
 * Created by mtutaj on 3/14/2018.
 */
public class AgrLoader {

    XdbIdDAO xdao = new XdbIdDAO();

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
