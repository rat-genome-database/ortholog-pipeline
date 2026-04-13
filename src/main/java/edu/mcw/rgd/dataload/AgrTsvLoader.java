package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by mtutaj on 3/14/2018.
 */
public class AgrTsvLoader {

    OrthologRelationDao dao = new OrthologRelationDao();

    Logger log = LogManager.getLogger("agrStatus");

    private String allianceFile;
    private String obsoleteOrthologsDeleteThreshold;
    private Set<String> processedSpecies;

    Map <String, List<Integer>> allianceIdToGeneRgdIdMap;

    boolean qcSymbolsForHumanGenes = true;

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

        int initialOrthologCount = dao.getAgrOrthologCount();
        log.info("initial ortholog count: "+initialOrthologCount);

        allianceIdToGeneRgdIdMap = dao.getAllianceIdToGeneRgdIdMap();
        log.info("loaded map of alliance-id to active-gene-in-rgd: "+allianceIdToGeneRgdIdMap.size());

        // list of all lines with species that are in RGD
        List<LineData> validLines = loadLines();

        CounterPool counters = new CounterPool();
        Set<Integer> accXdbKeys = run2(validLines, counters);

        wrapUp(time0, initialOrthologCount, accXdbKeys);

        log.info(counters.dumpAlphabetically());

        log.info("===== OK =====   elapsed "+Utils.formatElapsedTime(startTime, System.currentTimeMillis()));
    }

    List<LineData> loadLines() throws Exception {

        Set<Integer> processedSpeciesTypeKeys = new HashSet<>();
        for( String processedSpeciesName: getProcessedSpecies() ) {
            int sp = SpeciesType.parse(processedSpeciesName);
            if( sp>0 ) {
                processedSpeciesTypeKeys.add(sp);
            }
        }

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
                d.geneSymbol1 = parseSymbol(cols[1]);
                d.taxonId1 = cols[2]; // f.e. "NCBITaxon:9606"
                d.curie2 = cols[4];
                d.geneSymbol2 = parseSymbol(cols[5]);
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

        return validLines;
    }

    public Set<Integer> run2(List<LineData> list, CounterPool counters) {

        AtomicInteger inserted = new AtomicInteger(0);
        AtomicInteger updated = new AtomicInteger(0);
        AtomicInteger skipped = new AtomicInteger(0);

        Set<Integer> accXdbKeys = new HashSet<>();

        AtomicInteger remainingLines = new AtomicInteger(list.size());

        list.stream().parallel().forEach( d -> {

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

                        int r = dao.upsertAgrOrtholog(g1RgdId, g2RgdId, confidence, isBestScore, isBestRevScore, d.algorithms);
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
        });
        //}

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
        BufferedReader in = Utils.openReaderUtf8(localFile);
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
        fd.setAppendDateStamp(true);
        fd.setExternalFile(getAllianceFile());
        fd.setLocalFile("data/alliance.tsv");
        String localFile = fd.downloadNew();
        log.info("downloaded "+getAllianceFile());
        return localFile;
    }

    void wrapUp(Date time0, int initialOrthologCount, Set<Integer> curieAccXdbKeys) throws Exception {

        // delete stale orthologs
        int orthologCount = dao.getAgrOrthologCount();
        log.info("current ortholog count: "+orthologCount);

        int staleOrthologCount = dao.logStaleAgrOrthologs(time0);
        log.info(" stale ortholog count: "+staleOrthologCount);

        // convert delete-threshold-in-percent to a number
        int maxObsoleteOrthologsToBeDeleted;
        if( obsoleteOrthologsDeleteThreshold.endsWith("%") ) {
            int threshold = Integer.parseInt(obsoleteOrthologsDeleteThreshold.substring(0, obsoleteOrthologsDeleteThreshold.length()-1));
            maxObsoleteOrthologsToBeDeleted = (threshold * orthologCount) / 100;
        } else {
            maxObsoleteOrthologsToBeDeleted = Integer.parseInt(obsoleteOrthologsDeleteThreshold);
        }

        int newOrthologCount = orthologCount - staleOrthologCount;
        if( Math.abs(newOrthologCount-initialOrthologCount) > maxObsoleteOrthologsToBeDeleted ) {
            log.warn("*** WARN *** Cannot delete more than "+obsoleteOrthologsDeleteThreshold+" ("+maxObsoleteOrthologsToBeDeleted+") of orthologs!");
        } else {
            int staleRowsDeleted = dao.deleteStaleAgrOrthologs(time0);
            log.info("stale rows deleted from AGR_ORTHOLOGS: " + staleRowsDeleted);
        }

        orthologCount = dao.getAgrOrthologCount();
        log.info("final ortholog count: "+orthologCount);

        // QC AGR_GENE  curies
        dao.qcCuries(time0, curieAccXdbKeys, log);
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

                int issues = validateGeneSymbol(speciesTypeKey, geneSymbol, geneId, geneRgdId, counters);
                if( issues!=0 ) {
                    counters.add("*** GENE SYMBOL PROBLEMS", issues);
                }
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
            log.debug("trying to insert agr gene: species="+speciesTypeKey+", geneSymbol=["+geneSymbol+"], geneId="+geneId);
            Gene gene = dao.insertAgrGene(speciesTypeKey, geneSymbol, geneId);
            if( gene!=null ) {
                geneRgdId = gene.getRgdId();
                counters.increment("RESOLVE_GENE by insertion");
            }
        }
        return geneRgdId;
    }

    int validateGeneSymbol(int speciesTypeKey, String geneSymbol, String geneId, int geneRgdId, CounterPool counters) throws Exception {

        if( !this.qcSymbolsForHumanGenes ) {
            return 0;
        }

        int issueCount = 0;

        if( speciesTypeKey==SpeciesType.HUMAN ) {

            List<Gene> genes = dao.getGenesByXdbId(geneId, OrthologRelationDao.XDB_KEY_AGR_GENE);
            for( Gene g: genes ) {
                if( !Utils.stringsAreEqualIgnoreCase(g.getSymbol(), geneSymbol) ) {

                    if( !Utils.stringsAreEqualIgnoreCase(g.getEnsemblGeneSymbol(), geneSymbol) ) {
                        log.warn("WARNING! Alliance gene-id " + geneId + " points to gene RGD:" + g.getRgdId() + " [" + g.getSymbol() + "]");
                        log.warn("    Alliance had gene symbol [" + geneSymbol + "]");
                        issueCount++;
                    } else {
                        counters.increment("Alliance matches by Ensembl gene symbol");
                    }
                }
            }

            Gene g = dao.getGeneByRgdId(geneRgdId);
            if( g!=null ) {
                if( !Utils.stringsAreEqualIgnoreCase(g.getSymbol(), geneSymbol) ) {

                    if( !Utils.stringsAreEqualIgnoreCase(g.getEnsemblGeneSymbol(), geneSymbol) ) {
                        log.warn("WARNING! Matching gene RGD:" + g.getRgdId() + " [" + g.getSymbol() + "]");
                        log.warn("    Alliance had gene symbol [" + geneSymbol + "]");
                        issueCount++;
                    } else {
                        counters.increment("Alliance matches by Ensembl gene symbol");
                    }
                }
            }
        }

        return issueCount;
    }

    String parseSymbol( String s ) {

        String s2 = s;
        byte[] b2 = s.getBytes(StandardCharsets.UTF_8);
        if( b2.length != s.length() ) {

            StringBuffer out = new StringBuffer();
            for( int i=0; i<s.length(); i++ ) {
                char c = s.charAt(i);
                int codePoint = s.codePointAt(i);
                if( codePoint==916 ) {
                    out.append("DELTA");
                } else if( codePoint==945 ) {
                    out.append("alpha");
                } else if( codePoint==946 ) {
                    out.append("beta");
                } else if( codePoint==947 ) {
                    out.append("gamma");
                } else if( codePoint==948 ) {
                    out.append("delta");
                } else if( codePoint==949 ) {
                    out.append("epsilon");
                } else if( codePoint==950 ) {
                    out.append("zeta");
                } else if( codePoint==951 ) {
                    out.append("lambda");
                } else if( codePoint==952 ) {
                    out.append("theta");
                } else if( codePoint==953 ) {
                    out.append("iota");
                } else if( codePoint==954 ) {
                    out.append("kappa");
                } else if( codePoint==955 ) {
                    out.append("lambda");
                } else if( codePoint==956 ) {
                    out.append("mu");
                } else if( codePoint==963 ) {
                    out.append("sigma");
                } else if( codePoint > 127 ) {
                    log.error("#### unhandled char with codepoint "+codePoint);
                } else {
                    out.append(c);
                }
            }

            s2 = out.toString();
            //System.out.println(s+" "+s2);
        }
        return s2;
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
