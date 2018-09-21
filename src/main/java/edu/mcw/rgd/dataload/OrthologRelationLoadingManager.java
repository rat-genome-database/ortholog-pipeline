package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;
import edu.mcw.rgd.log.RGDSpringLogger;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;


/*
 * @author dli
 * @since Mar 23, 2007
 * <p>
 * load ortholog HGNC HCOP data from HGNC, for rat, mouse, human and dog
 * and load ortholog data from NCBI for every species
 */
public class OrthologRelationLoadingManager {

    OrthologRelationParser orthologRelationParser;
    OrthologRelationLoader loader;
    private String xrefDataSrc = "HGNC";

    long runTime = 0;
    int countOfRatGenesWithoutOrthologs = 0;
    RGDSpringLogger rgdLogger;

    protected final Logger status = Logger.getLogger("status");
    protected final Logger process = Logger.getLogger("process");
    protected final Logger logNoOrthologs = Logger.getLogger("NoOrthologs");
    protected final Logger logCrossLinked = Logger.getLogger("CrossLinkedOrthologs");

    static private OrthologRelationLoadingManager _instance;
    private String version;

    static public OrthologRelationLoadingManager getInstance() {
        return _instance;
    }

    public static void main(String[] args) {

        // parse cmd line params
        int speciesTypeKey = SpeciesType.ALL;
        boolean fixXrefDataSet = false;
        for( int i=0; i<args.length; i++ ) {
            switch(args[i]) {
                case "--fixXRefDataSet":
                    fixXrefDataSet = true;
                    break;
                case "--species":
                    speciesTypeKey = SpeciesType.parse(args[++i]);
                    break;
            }
        }
        if( speciesTypeKey==SpeciesType.ALL || speciesTypeKey==SpeciesType.HUMAN ) {
            System.out.println("ERROR: --cmdline parameter --species not specified or --species==human");
            System.exit(-1);
            return;
        }

        // run the instance
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        _instance = (OrthologRelationLoadingManager) (bf.getBean("orthologRelationLoadingManager"));
        try {
            if( fixXrefDataSet ) {
                _instance.fixXrefDataSet();
            } else {
                _instance.run(speciesTypeKey);
            }
        }
        catch( Exception e ) {
            _instance.process.error("CRITICAL ERROR", e);
            e.printStackTrace();
        }
    }
        
    public void run(int speciesTypeKey) throws Exception {

        String species = SpeciesType.getCommonName(speciesTypeKey).toUpperCase();

        status.info(getVersion());

        long startMilisec=System.currentTimeMillis();

        int oldOrthoCount = loader.dao.getOrthologCount(speciesTypeKey, SpeciesType.HUMAN);
        int oldAssocCount = loader.dao.getAssociationCountByType("weak_ortholog", speciesTypeKey, SpeciesType.HUMAN);
        status.info(species+"-HUMAN OLD ORTHOLOG COUNT = "+oldOrthoCount);
        status.info(species+"-HUMAN OLD ASSOCIATION COUNT = "+oldAssocCount);

        List<OrthologRelation> orthos = orthologRelationParser.parse(speciesTypeKey);
        loader.run(orthos, speciesTypeKey);

        long endMilisec=System.currentTimeMillis();
        runTime=(endMilisec - startMilisec)/1000;
        logRgd();

        int newOrthoCount = loader.dao.getOrthologCount(speciesTypeKey, SpeciesType.HUMAN);
        int newAssocCount = loader.dao.getAssociationCountByType("weak_ortholog", speciesTypeKey, SpeciesType.HUMAN);

        int diffOrthoCount = newOrthoCount - oldOrthoCount;
        int diffAssocCount = newAssocCount - oldAssocCount;
        NumberFormat plusMinusNF = new DecimalFormat(" +#; -#");

        String diffOrtho = diffOrthoCount!=0 ? "     difference: "+ plusMinusNF.format(diffOrthoCount) : "";
        String diffAssoc = diffAssocCount!=0 ? "     difference: "+ plusMinusNF.format(diffAssocCount) : "";

        status.info(species+"-HUMAN NEW ORTHOLOG COUNT = "+newOrthoCount + diffOrtho);
        status.info(species+"-HUMAN NEW ASSOCIATION COUNT = "+newAssocCount + diffAssoc);

        status.info("Processing time: "+Utils.formatElapsedTime(startMilisec, endMilisec));
        status.info("===");
    }
   
    public void logRgd() throws Exception {
        process.info("Logging into RGD database!");
        int unmatched = loader.getUnMatchedC() + orthologRelationParser.getSkippedRelations() ;
        
        rgdLogger.log("orthologRelations","totalNumberRelationsProcessed", orthologRelationParser.getTotalRelations());
        rgdLogger.log("orthologRelations","relationsUnmatched", unmatched);
        rgdLogger.log("orthologRelations","relationsMatchedWithdrawnGene", loader.getWithdrawnC());
        rgdLogger.log("orthologRelations","relationsInserted", loader.getRowsInserted());
        rgdLogger.log("orthologRelations","relationsDeleted", loader.getRowsDeleted());
        rgdLogger.log("orthologRelations","relationsMultipleMatched", loader.getMultipleMatchC());
        if( countOfRatGenesWithoutOrthologs!=0 ) {
            rgdLogger.log("orthologRelations","ratGenesWithoutOrthologs", countOfRatGenesWithoutOrthologs);
        }
        rgdLogger.log("orthologRelations","timeToExecute", getRunTime());
    }

    void dumpCrossLinkedOrthologs() throws Exception {

        logCrossLinked.warn("=====");
        logCrossLinked.warn("Report generated only for active rat genes");
        logCrossLinked.warn("=====");
        logCrossLinked.warn("RGD_ID\tGENE_SYMBOL");
        for( String crossLinkedInfo: loader.getOrthologRelationDao().getCrossLinkedOrthologs(SpeciesType.RAT) ) {
            logCrossLinked.warn(crossLinkedInfo);
        }
    }

    /**
     *
     * @return count of active rat genes without orthologs
     * @throws Exception
     */
    int dumpGenesWithoutOrthologs() throws Exception {

        logNoOrthologs.info("gene type> space separated list of GENE_SYMBOL:RGD_ID");

        List<Gene> genes = loader.getOrthologRelationDao().getGenesWithoutOrthologs(SpeciesType.RAT);
        // sort genes by symbol
        Collections.sort(genes, new Comparator<Gene>() {
            public int compare(Gene o1, Gene o2) {
                return Utils.stringsCompareToIgnoreCase(o1.getSymbol(), o2.getSymbol());
            }
        });

        int count = genes.size();
        while( !genes.isEmpty() ) {
            Gene gene = genes.get(0);
            String geneTypeLc = Utils.defaultString(gene.getType());

            // dump genes without orthologs for genes of given type
            String geneList = dumpGenesWithoutOrthologs(geneTypeLc, genes);
            if( !geneTypeLc.equals("allele") && !geneTypeLc.equals("splice") ) {
                logNoOrthologs.info(geneTypeLc + "> " + geneList);
            }
        }
        return count;
    }

    String dumpGenesWithoutOrthologs(String geneType, List<Gene> genes) {

        StringBuilder buf = new StringBuilder();
        Iterator<Gene> it = genes.iterator();
        while( it.hasNext() ) {
            Gene gene = it.next();
            if( !Utils.stringsAreEqual(geneType, gene.getType()) )
                continue;

            buf.append(" ").append(gene.getSymbol()).append(":").append(gene.getRgdId());
            it.remove();
        }
        return buf.toString();
    }

    void fixXrefDataSet() throws Exception {

        OrthologRelationDao dao = new OrthologRelationDao();

        // fix XREF_DATA_SET in orthologs
        int rowsAffected = dao.fixXrefDataSetInOrthologs();
        status.info("Fixed XREF_DATA_SET for orthologs: "+rowsAffected);

        rowsAffected = dao.fixXrefDataSetInAssociations();
        status.info("Fixed XREF_DATA_SET for associations: "+rowsAffected);
    }

    public OrthologRelationLoader getOrthologRelationLoader() {
        return loader;
    }
    
    public void setOrthologRelationLoader( OrthologRelationLoader loader) {
        this.loader = loader;
    }
    
    public OrthologRelationParser getOrthologRelationParser() {
        return orthologRelationParser;
    }
    public void setOrthologRelationParser(OrthologRelationParser orthologRelationParser) {
        this.orthologRelationParser = orthologRelationParser;
    }
    
    public long getRunTime() {
        return runTime;
    }    
    
    public RGDSpringLogger getRgdLogger() {
        return rgdLogger;
    }
    public void setRgdLogger(RGDSpringLogger rgdLogger) {
        this.rgdLogger = rgdLogger;
    }

    public void setXrefDataSrc(String xrefDataSrc) {
        this.xrefDataSrc = xrefDataSrc;
    }

    public String getXrefDataSrc() {
        return xrefDataSrc;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }
}
