package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.CounterPool;
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
    private Set<String> processedSpecies;

    Map <String, List<Integer>> allianceIdToGeneRgdIdMap;

    public void run() throws Exception {

        try {
            run2();
        } catch( Exception e ) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    public void run2() throws Exception {

        long startTime = System.currentTimeMillis();

        // note: there could be disparity between oracle server time and the time on machine the pipeline is running,
        //  so to prevent deletion of valid orthologs (classified by the code as 'stale orthologs'),
        //  we use the as the cutoff timestamp the machine timestamp minus one hour
        Date time0 = Utils.addHoursToDate(new Date(), -1); // time0 = current day-and-time less 1 hour

        Set<Integer> processedSpeciesTypeKeys = new HashSet<>();
        for( String processedSpeciesName: getProcessedSpecies() ) {
            int sp = SpeciesType.parse(processedSpeciesName);
            if( sp>0 ) {
                processedSpeciesTypeKeys.add(sp);
            }
        }

        int initialOrthologCount = getOrthologCount();
        log.info("initial ortholog count: "+initialOrthologCount);

        allianceIdToGeneRgdIdMap = dao.getAllianceIdToGeneRgdIdMap();
        log.info("loaded map of alliance-id to active-gene-in-rgd: "+allianceIdToGeneRgdIdMap.size());

        AtomicInteger skipped = new AtomicInteger(0);

        // list of all lines with species that are in RGD
        List<LineData> validLines = new ArrayList<>();
        {
            List<String> lines = loadLinesFromTsvFile();
            Collections.shuffle(lines);
            log.info("  lines read from tsv file: " + lines.size());

            Set<String> agrCuries = new HashSet<>();

            for( String line: lines ) {

                // data line, f.e.
                // HGNC:25018	TMEM216	NCBITaxon:9606	Homo sapiens	FB:FBgn0037615	CG11760	NCBITaxon:7227	Drosophila melanogaster	Ensembl Compara|InParanoid|Roundup|PhylomeDB|OMA|PANTHER|OrthoFinder|OrthoInspector|TreeFam|Hieranoid	10	10	Yes	Yes
                String[] cols = line.split("[\\t]", -1);
                String taxonName1 = cols[3]; // f.e. "Homo sapiens"
                String taxonName2 = cols[7]; // f.e. "Homo sapiens"

                int speciesTypeKey1 = SpeciesType.parse(taxonName1);
                int speciesTypeKey2 = SpeciesType.parse(taxonName2);

                if( !processedSpeciesTypeKeys.contains(speciesTypeKey1) ||
                        !processedSpeciesTypeKeys.contains(speciesTypeKey2) ) {

                    skipped.incrementAndGet();
                    continue;
                }

                LineData d = new LineData();
                d.taxonName1 = taxonName1;
                d.taxonName2 = taxonName2;
                d.speciesTypeKey1 = speciesTypeKey1;
                d.speciesTypeKey2 = speciesTypeKey2;

                d.curie1 = cols[0];
                d.geneSymbol1 = cols[1];
                d.taxonId1 = cols[2]; // f.e. "NCBITaxon:9606"
                d.curie2 = cols[4];
                d.geneSymbol2 = cols[5];
                d.taxonId2 = cols[6]; // f.e. "NCBITaxon:9606"

                d.algorithms = sortAlgorithmsStr(cols[8]);
                d.algorithmsMatch = cols[9]; // f.e. "6"
                d.outOfAlgorithms = cols[10]; // f.e. "9"
                d.isBestScoreStr = cols[11]; // True of False
                d.isBestRevScoreStr = cols[12]; // True or False
                validLines.add(d);

                agrCuries.add(d.curie1);
                agrCuries.add(d.curie2);
            }

            log.info("  incoming skipped  orthologs: "+skipped+" (species skipped from processing)");
            log.info("  incoming unique agr curies: "+agrCuries.size());
        }

        CounterPool counters = new CounterPool();
        Set<Integer> accXdbKeys = run2(validLines, counters);

        wrapUp(time0, initialOrthologCount, accXdbKeys);

        log.info(counters.dumpAlphabetically());

        log.info("===== OK =====   elapsed "+Utils.formatElapsedTime(startTime, System.currentTimeMillis()));
    }

    public Set<Integer> run2(List<LineData> list, CounterPool counters) {

        AtomicInteger inserted = new AtomicInteger(0);
        AtomicInteger updated = new AtomicInteger(0);
        AtomicInteger skipped = new AtomicInteger(0);

        Set<Integer> accXdbKeys = new HashSet<>();

        AtomicInteger remainingLines = new AtomicInteger(list.size());

        //for( LineData d: list ) {
        list.parallelStream().forEach( d -> {

            try {

                boolean wasProcessedWithoutExceptions = false;
                while (!wasProcessedWithoutExceptions) {
                    try {
                        int g1RgdId = resolveGene(d.speciesTypeKey1, d.geneSymbol1, d.curie1, counters);
                        if (g1RgdId == 0) {
                            log.warn("WARN: cannot resolve gene [" + d.geneSymbol1 + "] [" + d.curie1 + "] " + d.taxonName1);
                            skipped.incrementAndGet();
                            break;
                        }
                        resolveCurie(accXdbKeys, d.curie1, g1RgdId);

                        int g2RgdId = resolveGene(d.speciesTypeKey2, d.geneSymbol2, d.curie2, counters);
                        if (g2RgdId == 0) {
                            log.warn("WARN: cannot resolve gene [" + d.geneSymbol2 + "] [" + d.curie2 + "] " + d.taxonName2);
                            skipped.incrementAndGet();
                            break;
                        }
                        resolveCurie(accXdbKeys, d.curie2, g2RgdId);

                        boolean isBestScore = d.isBestScoreStr.equals("Yes");
                        boolean isBestRevScore = d.isBestRevScoreStr.equals("Yes");

                        String confidence = "stringent";

                        int r = updateDb(g1RgdId, g2RgdId, confidence, isBestScore, isBestRevScore, d.algorithms);
                        if (r == 1) {
                            inserted.incrementAndGet();
                        } else {
                            updated.incrementAndGet();
                        }

                        wasProcessedWithoutExceptions = true;
                        remainingLines.decrementAndGet();

                    } catch (org.springframework.dao.DuplicateKeyException ignore) {
                        counters.increment("org.springframework.dao.DuplicateKeyException");
                        log.debug("dao.DuplicateKeyException " + d.curie1 + " " + d.curie2 + "    remaining lines: " + remainingLines.get());
                    }
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            //});
        }

        log.info("inserted orthologs: "+inserted);
        log.info("updated  orthologs: "+updated);
        log.info("skipped  orthologs: "+skipped);

        return accXdbKeys;
    }

    void resolveCurie(Set<Integer> xdbIdKeys, String curie, int rgdId) throws Exception {

        int xdbKey = 63; // AGR_GENE
        List<XdbId> xdbIds = dao.getXdbIdsByRgdId(xdbKey, rgdId);
        for( XdbId xdbId: xdbIds ) {
            if( xdbId.getAccId().equals(curie) ) {
                synchronized (xdbIdKeys) {
                    xdbIdKeys.add(xdbId.getKey());
                }
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
    }

    int getOrthologCount() throws Exception {
        String sql = "SELECT COUNT(0) FROM agr_orthologs";
        return xdao.getCount(sql);
    }

    synchronized int resolveGene(int speciesTypeKey, String geneSymbol, String geneId, CounterPool counters) throws Exception {

        int geneRgdId = 0;

        /* original code -- slow
            // first resolve by xdb id
            int xdbKey = 63; // AGR_GENE
            List<Gene> genes = dao.getGenesByXdbId(geneId, xdbKey);
            if( genes.size()>1 ) {
                Collections.sort(genes, new Comparator<Gene>() {
                    @Override
                    public int compare(Gene o1, Gene o2) {
                        return o1.getRgdId() - o2.getRgdId();
                    }
                });
                log.warn("  *** MULTIPLE GENES FOR "+geneId+"  RGD:"+Utils.concatenate(",RGD:", genes, "getRgdId")+";  picking the first one");
            }
            if( !genes.isEmpty() ) {
                geneRgdId = genes.get(0).getRgdId();
            }
        */
        List<Integer> geneRgdIds = this.allianceIdToGeneRgdIdMap.get(geneId);
        if( geneRgdIds!=null ) {
            if (geneRgdIds.size() > 1) {
                log.warn("  *** MULTIPLE GENES FOR " + geneId + "  " + Utils.concatenate(geneRgdIds, ",") + ";  picking the first one");
            }
            if (!geneRgdIds.isEmpty()) {
                geneRgdId = geneRgdIds.get(0);
                counters.increment("RESOLVE_GENE from map");
            }
        }

        // second: special resolution for rat,mouse,human
        if( geneRgdId==0 ) {
            if( speciesTypeKey==SpeciesType.RAT ) {
                int geneRgdIdStr = Integer.parseInt(geneId.substring(4));
                Gene gene = dao.getGeneByRgdId(geneRgdIdStr);
                if( gene!=null ) {
                    geneRgdId = gene.getRgdId();
                }
            }
            else if( speciesTypeKey==SpeciesType.MOUSE ) {
                List<Gene> genes = dao.getGenesByXdbId(geneId, XdbId.XDB_KEY_MGD);
                if( !genes.isEmpty() ) {
                    geneRgdId = genes.get(0).getRgdId();
                }
            }
            else if( speciesTypeKey==SpeciesType.HUMAN ) {
                List<Gene> genes = dao.getGenesByXdbId(geneId, XdbId.XDB_KEY_HGNC);
                if( !genes.isEmpty() ) {
                    geneRgdId = genes.get(0).getRgdId();
                }
            }

            if( geneRgdId!=0 ) {
                // gene resolved: insert AGR_GENE xdbid
                dao.insertAgrGeneXdbId(geneRgdId, geneId);
                counters.increment("RESOLVE_GENE by rat/mouse/human special id");
            }
        }

        // second, resolve by gene symbol
        if( geneRgdId==0 ) {
            Gene gene = dao.getGeneBySymbol(geneSymbol, speciesTypeKey, log);
            if( gene!=null ) {
                geneRgdId = gene.getRgdId();
                // gene resolved: insert AGR_GENE xdbid
                dao.insertAgrGeneXdbId(geneRgdId, geneId);
                counters.increment("RESOLVE_GENE by gene symbol");
            }
        }

        // if gene still not found, insert it
        if( geneRgdId==0 ) {
            Gene gene = dao.insertAgrGene(speciesTypeKey, geneSymbol, geneId);
            if( gene!=null ) {
                geneRgdId = gene.getRgdId();
                counters.increment("RESOLVE_GENE by insertion");
            }
        }
        return geneRgdId;
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

    public void setProcessedSpecies(Set<String> processedSpecies) {
        this.processedSpecies = processedSpecies;
    }

    public Set<String> getProcessedSpecies() {
        return processedSpecies;
    }

    class LineData {

        String curie1;
        String geneSymbol1;
        String taxonId1; // f.e. "NCBITaxon:9606"
        String taxonName1; // f.e. "Homo sapiens"
        String curie2;
        String geneSymbol2;
        String taxonId2; // f.e. "NCBITaxon:9606"
        String taxonName2; // f.e. "Homo sapiens"

        String algorithms;
        String algorithmsMatch; // f.e. "6"
        String outOfAlgorithms; // f.e. "9"
        String isBestScoreStr; // True of False
        String isBestRevScoreStr; // True or False

        // derived fields
        int speciesTypeKey1;
        int speciesTypeKey2;
    }
}
