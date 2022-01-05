package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.object.MappingSqlQuery;

import java.io.BufferedReader;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by mtutaj on 3/14/2018.
 */
public class AgrTsvLoader {

    XdbIdDAO xdao = new XdbIdDAO();
    OrthologRelationDao dao = new OrthologRelationDao();

    Logger log = LogManager.getLogger("agrStatus");
    Logger logDel = LogManager.getLogger("deletedAgrOrthologs");
    Logger logIns = LogManager.getLogger("insertedAgrOrthologs");

    private String allianceFile;
    private String obsoleteOrthologsDeleteThreshold;

    public void run() throws Exception {

        try {
            run2();
        } catch( Exception e ) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    public void run2() throws Exception {

        // note: there could be disparity between oracle server time and the time on machine the pipeline is running,
        //  so to prevent deletion of valid orthologs (classified by the code as 'stale orthologs'),
        //  we use the as the cutoff timestamp the machine timestamp minus one hour
        Date time0 = Utils.addHoursToDate(new Date(), -1); // time0 = current day-and-time less 1 hour

        int initialOrthologCount = getOrthologCount();
        log.info("initial ortholog count: "+initialOrthologCount);

        AtomicInteger inserted = new AtomicInteger(0);
        AtomicInteger updated = new AtomicInteger(0);

        List<String> lines = loadLinesFromTsvFile();
        Collections.shuffle(lines);
        log.info("  LINES SHUFFLED");

        Set<Integer> accXdbKeys = new HashSet<>();

        //int i = 0;
        for( String line: lines ) {

            //log.debug((++i)+".");

            // data line, f.e.
            // HGNC:25018	TMEM216	NCBITaxon:9606	Homo sapiens	FB:FBgn0037615	CG11760	NCBITaxon:7227	Drosophila melanogaster	Ensembl Compara|InParanoid|Roundup|PhylomeDB|OMA|PANTHER|OrthoFinder|OrthoInspector|TreeFam|Hieranoid	10	10	Yes	Yes
            String[] cols = line.split("[\\t]", -1);

            String curie1 = cols[0];
            String geneSymbol1 = cols[1];
            String taxonId1 = cols[2]; // f.e. "NCBITaxon:9606"
            String taxonName1 = cols[3]; // f.e. "Homo sapiens"
            String curie2 = cols[4];
            String geneSymbol2 = cols[5];
            String taxonId2 = cols[6]; // f.e. "NCBITaxon:9606"
            String taxonName2 = cols[7]; // f.e. "Homo sapiens"

            String algorithms = sortAlgorithmsStr(cols[8]);
            String algorithmsMatch = cols[9]; // f.e. "6"
            String outOfAlgorithms = cols[10]; // f.e. "9"
            String isBestScoreStr = cols[11]; // True of False
            String isBestRevScoreStr = cols[12]; // True or False

            int speciesTypeKey1 = SpeciesType.parse(taxonName1);
            Gene g1 = resolveGene(speciesTypeKey1, geneSymbol1, curie1);
            if( g1==null ) {
                log.warn("WARN: cannot resolve gene ["+geneSymbol1+"] ["+curie1+"] "+taxonName1);
                continue;
            }
            resolveCurie(accXdbKeys, curie1, g1.getRgdId());

            int speciesTypeKey2 = SpeciesType.parse(taxonName2);
            Gene g2 = resolveGene(speciesTypeKey2, geneSymbol2, curie2);
            if( g2==null ) {
                log.warn("WARN: cannot resolve gene ["+geneSymbol2+"] ["+curie2+"] "+taxonName2);
                continue;
            }
            resolveCurie(accXdbKeys, curie2, g2.getRgdId());

            boolean isBestScore = isBestScoreStr.equals("Yes");
            boolean isBestRevScore = isBestRevScoreStr.equals("Yes");

            String confidence = "stringent";

            int r = updateDb(g1.getRgdId(), g2.getRgdId(), confidence, isBestScore, isBestRevScore, algorithms);
            if( r==1 ) {
                inserted.incrementAndGet();
            } else {
                updated.incrementAndGet();
            }
        }

        log.info("inserted orthologs: "+inserted);
        log.info("updated  orthologs: "+updated);

        wrapUp(time0, initialOrthologCount, accXdbKeys);
    }

    void resolveCurie(Set<Integer> xdbIdKeys, String curie, int rgdId) throws Exception {

        int xdbKey = 63; // AGR_GENE
        List<XdbId> xdbIds = dao.getXdbIdsByRgdId(xdbKey, rgdId);
        for( XdbId xdbId: xdbIds ) {
            if( xdbId.getAccId().equals(curie) ) {
                xdbIdKeys.add(xdbId.getKey());
                return;
            }
        }

        //throw new Exception("resolveCurie for XDB_KEY=63 failed: "+curie+", RGD:"+rgdId);
    }

    /// sort '|'-separated algorithm string, to present the data in nice way
    String sortAlgorithmsStr(String s) {
        String[] arr = s.split("\\|");
        Set<String> sorted = new TreeSet<>(Arrays.asList(arr));
        return Utils.concatenate(sorted, "|");
    }

    List<String> loadLinesFromTsvFile() throws Exception {

        String localFile = downloadTsvFile();
        BufferedReader in = Utils.openReader(localFile);
        int dataLinesRead = 0;
        List<String> lines = new ArrayList<>();

        String line;
        while( (line=in.readLine()) != null ) {

            // skip comment lines
            if( line.startsWith("#") ) {
                continue;
            }

            // expected header line:
            // Gene1ID	Gene1Symbol	Gene1SpeciesTaxonID	Gene1SpeciesName	Gene2ID	Gene2Symbol	Gene2SpeciesTaxonID	Gene2SpeciesName	Algorithms	AlgorithmsMatch	OutOfAlgorithms	IsBestScore	IsBestRevScore
            if( dataLinesRead++ == 0 ) {
                // skip header line
                continue;
            }

            lines.add(line);
        }
        in.close();

        log.info("incoming data lines: "+dataLinesRead);

        return lines;
    }

    String downloadTsvFile() throws Exception {
        FileDownloader fd = new FileDownloader();
        fd.setUseCompression(true);
        fd.setPrependDateStamp(true);
        fd.setExternalFile(getAllianceFile());
        fd.setLocalFile("data/ORTHOLOGY-ALLIANCE_COMBINED.tsv");
        String localFile = fd.downloadNew();
        log.info("downloaded "+getAllianceFile());
        return localFile;
    }

    void wrapUp(Date time0, int initialOrthologCount, Set<Integer> curieAccXdbKeys) throws Exception {

        // delete stale orthologs
        int orthologCount = getOrthologCount();
        log.info("current ortholog count: "+orthologCount);

        String sql = "SELECT * FROM agr_orthologs WHERE last_update_date<?";
        MappingSqlQuery q = new MappingSqlQuery(xdao.getDataSource(), sql) {
            @Override
            protected Object mapRow(ResultSet rs, int rowNum) throws SQLException {
                int geneRgdId1 = rs.getInt("gene_rgd_id_1");
                int geneRgdId2 = rs.getInt("gene_rgd_id_2");
                String conf = rs.getString("confidence");
                String isBest = rs.getString("is_best_score");
                String isBestRev = rs.getString("is_best_rev_score");
                String methodsMatched = rs.getString("methods_matched");
                Date createdDate = rs.getTimestamp("created_date");
                Date lastUpdateDate = rs.getTimestamp("last_update_date");

                logDel.debug("RGD1:"+geneRgdId1+"|RGD2:"+geneRgdId2+"|"+conf+
                    "|IS_BEST:"+isBest+"|IS_BEST_REV:"+isBestRev+"|METHODS_MATCHED:"+methodsMatched+
                    "|CREATED:"+createdDate+"|LAST_UPDATED:"+lastUpdateDate);
                return null;
            }
        };
        q.declareParameter(new SqlParameter(Types.TIMESTAMP));
        q.compile();
        List staleOrthologs = q.execute(time0);
        log.info(" stale ortholog count: "+staleOrthologs.size());

        // convert delete-threshold-in-percent to a number
        int maxObsoleteOrthologsToBeDeleted;
        if( obsoleteOrthologsDeleteThreshold.endsWith("%") ) {
            int threshold = Integer.parseInt(obsoleteOrthologsDeleteThreshold.substring(0, obsoleteOrthologsDeleteThreshold.length()-1));
            maxObsoleteOrthologsToBeDeleted = (threshold * orthologCount) / 100;
        } else {
            maxObsoleteOrthologsToBeDeleted = Integer.parseInt(obsoleteOrthologsDeleteThreshold);
        }

        int newOrthologCount = orthologCount - staleOrthologs.size();
        if( Math.abs(newOrthologCount-initialOrthologCount) > maxObsoleteOrthologsToBeDeleted ) {
            log.warn("*** WARN *** Cannot delete more than "+obsoleteOrthologsDeleteThreshold+" ("+maxObsoleteOrthologsToBeDeleted+") of orthologs!");
        } else {
            sql = "DELETE FROM agr_orthologs WHERE last_update_date<?";
            int staleRowsDeleted = xdao.update(sql, time0);
            log.info("stale rows deleted from AGR_ORTHOLOGS: " + staleRowsDeleted);
        }

        orthologCount = getOrthologCount();
        log.info("final ortholog count: "+orthologCount);

        // QC AGR_GENE  curies
        dao.qcCuries(time0, curieAccXdbKeys, log);

        log.info("===== OK =====   elapsed "+Utils.formatElapsedTime(time0.getTime(), System.currentTimeMillis()));
    }

    int getOrthologCount() throws Exception {
        String sql = "SELECT COUNT(0) FROM agr_orthologs";
        return xdao.getCount(sql);
    }

    Gene resolveGene(int speciesTypeKey, String geneSymbol, String geneId) throws Exception {

        Gene gene = null;

        // first resolve by xdb id
        int xdbKey = 63; // AGR_GENE
        List<Gene> genes = dao.getGenesByXdbId(geneId, xdbKey);
        if( genes.size()>1 ) {
            log.warn("  *** MULTIPLE GENES FOR "+geneId+"  RGD:"+Utils.concatenate(",RGD:", genes, "getRgdId")+";  picking the first one");
            gene = genes.get(0);
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
            gene = dao.getGeneBySymbol(geneSymbol, speciesTypeKey, log);
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

    int updateDb(int rgdId1, int rgdId2, String confidence, boolean isBestScore, boolean isBestRevScore, String methodsMatched) throws Exception {

        String bestScore = isBestScore?"Y":"N";
        String bestRevScore = isBestRevScore?"Y":"N";

        String sql1 = "SELECT COUNT(0) FROM agr_orthologs WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";
        int dataExists = xdao.getCount(sql1, rgdId1, rgdId2, methodsMatched);
        if( dataExists==0 ) {
            // insert
            final String sql = "INSERT INTO agr_orthologs (gene_rgd_id_1, gene_rgd_id_2, confidence, is_best_score, "+
                    "is_best_rev_score, methods_matched, last_update_date) "+
                    "VALUES(?,?,?,?,?,?,SYSDATE)";
            xdao.update(sql, rgdId1, rgdId2, confidence, bestScore, bestRevScore, methodsMatched);

            logIns.debug("RGD1:"+rgdId1+"|RGD2:"+rgdId2+"|"+confidence+
                    "|IS_BEST:"+bestScore+"|IS_BEST_REV:"+bestRevScore+"|METHODS_MATCHED:"+methodsMatched);
            return 1;
        } else{
            final String sql = "UPDATE agr_orthologs SET confidence=?, is_best_score=?, is_best_rev_score=?, " +
                    "last_update_date=SYSDATE " +
                    "WHERE gene_rgd_id_1=? AND gene_rgd_id_2=? AND methods_matched=?";

            xdao.update(sql, confidence, bestScore, bestRevScore, rgdId1, rgdId2, methodsMatched);
            return 0;
        }
    }

    public void setAllianceFile(String allianceFile) {
        this.allianceFile = allianceFile;
    }

    public String getAllianceFile() {
        return allianceFile;
    }

    public void setObsoleteOrthologsDeleteThreshold(String obsoleteOrthologsDeleteThreshold) {
        this.obsoleteOrthologsDeleteThreshold = obsoleteOrthologsDeleteThreshold.trim();
    }
}
