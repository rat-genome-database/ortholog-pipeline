package edu.mcw.rgd.dataload;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.*;
import java.util.Map;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.RowMapper;

import javax.sql.DataSource;

/*
 * Created on Mar 28, 2007
 */
public class OrthologRelationDao {

    protected final Logger logger = Logger.getLogger("process");
    protected final Logger logInsertedOrthologs = Logger.getLogger("InsertedOrthologs");
    public final Logger logDeletedOrthologs = Logger.getLogger("DeletedOrthologs");

    private OrthologDAO orthologDAO = new OrthologDAO();
    private RGDManagementDAO rgdIdDAO = new RGDManagementDAO();
    private AssociationDAO associationDAO = new AssociationDAO();
    private XdbIdDAO xdbIdDAO = new XdbIdDAO();
    private GeneDAO geneDAO = associationDAO.getGeneDAO();

    private int directOrthologTypeKey;
    private int transitiveOrthologTypeKey;

    public OrthologRelationDao() {
        logger.info(orthologDAO.getConnectionInfo());
    }



    /**
     * get manual orthologs (with data source of 'RGD') for given source and dest species
     * @param srcRgdId rgd id
     * @param destSpeciesTypeKey species type key for destination rgd id
     * @return List of Ortholog objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<Ortholog> getManualOrthologs(int srcRgdId, int destSpeciesTypeKey) throws Exception {

        // load all orthologs for given src rgd id
        List<Ortholog> orthos = orthologDAO.getOrthologsForSourceRgdId(srcRgdId);
        // remove all non-manual orthologs
        Iterator<Ortholog> it = orthos.listIterator();
        while( it.hasNext() ) {
            Ortholog o = it.next();
            if( !Utils.stringsAreEqual(o.getXrefDataSrc(), "RGD") || o.getDestSpeciesTypeKey()!=destSpeciesTypeKey) {
                it.remove();
            }
        }
        return orthos;
    }

    public boolean deleteStaleOrtholog(Ortholog ortholog) throws Exception {

        // REQUIREMENT 1: never delete a manual ortholog
        if( Utils.stringsAreEqual(ortholog.getXrefDataSrc(), "RGD") ) {
            return false;
        }

        int srcRgdId = ortholog.getSrcRgdId();
        List<Ortholog> slist = orthologDAO.getOrthologsForSourceRgdId(srcRgdId);

        // remove from slist orthologs of different dest species type key
        Iterator<Ortholog> it = slist.iterator();
        while( it.hasNext() ) {
            Ortholog o = it.next();
            if( o.getDestSpeciesTypeKey()!=ortholog.getDestSpeciesTypeKey() )
                it.remove();
        }

        // slist contains now all orthologs with same rgd id and dest species type key
        // REQUIREMENT 2: cannot delete the only ortholog from database!!!
        if( slist.size()>1 ) {
            List<Ortholog> list = new ArrayList<>(1);
            list.add(ortholog);
            deleteOrthologs(list);
            return true;
        }
        return false;
    }

    /**
     * check if two genes are orthologous
     * @return lowest value of genetogene_key if there is an orthology between two genes; 0 otherwise
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int getKeyForMatchingOrtholog(Ortholog ortholog, List<Ortholog> deleteList) throws Exception {
        // check if these genes have multiple orthologs

        int srcRgdId = ortholog.getSrcRgdId();
        List<Ortholog> slist = orthologDAO.getOrthologsForSourceRgdId(srcRgdId);

        // remove from slist orthologs of different dest species type key
        Iterator<Ortholog> it = slist.iterator();
        while( it.hasNext() ) {
            Ortholog o = it.next();
            if( o.getDestSpeciesTypeKey()!=ortholog.getDestSpeciesTypeKey() )
                it.remove();
        }

        Ortholog existingOrtholog = null;
        if( slist.size()>1) {

            // conflict: we have multiple orthologs for srcRgdId and destSpeciesTypeKey in database already
            // sort those orthologs by priority, lowest priority first, and keep only one
            Collections.sort(slist, new Comparator<Ortholog>() {
                public int compare(Ortholog o1, Ortholog o2) {
                    return compareOrthologs(o1, o2);
                }
            });
            existingOrtholog = slist.remove(0);
            deleteList.addAll(slist);
        }
        else if( slist.size()==1 ) {
            existingOrtholog = slist.get(0);
        }

        if( existingOrtholog==null )
            return 0;
        else if( existingOrtholog.getDestRgdId()==ortholog.getDestRgdId() )
            return existingOrtholog.getKey();

        // there is no matching ortholog but
        // there is another ortholog in RGD for same srcRgdId and destSpeciesTypeKey;
        // see which ortholog is better match: the existing one in RGD or the being checked for
        int r = compareOrthologs(existingOrtholog, ortholog);
        if( r>0 ) {
            // existing in RGD ortholog is weaker than incoming ortholog: drop it from RGD
            List<Ortholog> orthologsForDel = new ArrayList<>(1);
            orthologsForDel.add(existingOrtholog);
            deleteOrthologs(orthologsForDel);

            // and add the incoming ortholog
            return 0;
        } else {
            // incoming ortholog is weaker than existing in-RGD ortholog; degrade incoming ortholog to association
            return -existingOrtholog.getKey();
        }
    }

    // methodCounters - array of how many times given best-fit method was used
    int compareOrthologs(Ortholog o1, Ortholog o2) {

        // pick ortholog with longest evidence
        int ev1count = getEvidenceCount(o1);
        int ev2count = getEvidenceCount(o2);
        if( ev1count != ev2count )
            return ev2count - ev1count;

        // pick relation with matching symbol
        String sourceGeneSymbol = getGeneSymbol(o1.getSrcRgdId());
        String dest1GeneSymbol = getGeneSymbol(o1.getDestRgdId());
        String dest2GeneSymbol = getGeneSymbol(o2.getDestRgdId());
        if( Utils.stringsCompareToIgnoreCase(sourceGeneSymbol, dest1GeneSymbol)==0 )
            return -1;
        if( Utils.stringsCompareToIgnoreCase(sourceGeneSymbol, dest2GeneSymbol)==0 )
            return 1;

        return dest2GeneSymbol.compareToIgnoreCase(dest1GeneSymbol);
    }

    int getEvidenceCount(Ortholog o) {

        int count = 0;
        String evidences = o.getXrefDataSet();
        if( evidences!=null ) {
            count = 1;
            for( int pos=0; (pos=evidences.indexOf(',', pos))>0; pos++ ) {
                count++;
            }
        }
        return count;
    }

    /**
     * check if two genes are orthologous
     * @param srcRgdId rgd id if first gene
     * @param destRgdId rgd id if second gene
     * @return lowest value of genetogene_key if there is an orthology between two genes; 0 otherwise
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int areGenesOrthologous(int srcRgdId, int destRgdId) throws Exception {
        return orthologDAO.areGenesOrthologous(srcRgdId, destRgdId);
    }

    /**
     * check if there is weak ortholog (association) between two genes
     * @param srcRgdId rgd id if first gene
     * @param destRgdId rgd id if second gene
     * @return true if there is a weak orthology between two genes; false otherwise
     * @throws Exception when unexpected error in spring framework occurs
     */
    public boolean areGenesAssociated(int srcRgdId, int destRgdId) throws Exception {
        List<Association> list = associationDAO.getAssociationsForMasterRgdId(srcRgdId, "weak_ortholog");
        for( Association assoc: list ) {
            if( assoc.getDetailRgdId()==destRgdId )
                return true;
        }
        return false;
    }

    public void checkComplementOrthologs(int speciesType1, int speciesType2, List<Association> forDeleteAssociations) throws Exception {

        int nondowngradedOrthologs = 0;
        int downgradedOrthologs = 0;
        int reverseOrthologMatch = 0;
        int reverseAssocMatch = 0;
        int reverseMissing = 0;
        int detachedAssociations = 0;

        // if there is an ortholog without matching reverse ortholog
        // then downgrade it into an association
        List<Ortholog> orthologsForDelete = new ArrayList<>();

        List<Ortholog> orthologs = orthologDAO.getAllOrthologs(speciesType1, speciesType2);

        for( Ortholog o: orthologs ) {
            // check reverse ortholog
            if( areGenesOrthologous(o.getDestRgdId(), o.getSrcRgdId())!=0 ) {
                reverseOrthologMatch++;
                continue; // reverse ortholog does exist
            }

            // there is no reverse ortholog! drop ortholog
            if( o.getXrefDataSrc().equals("RGD") ) {
                logger.error("attempt to delete manual ortholog without reverse ortholog "+o.dump("|"));
                nondowngradedOrthologs++;
            } else {
                orthologsForDelete.add(o);

                // and create association (weak ortholog) instead
                Association assoc = new Association();
                assoc.setAssocSubType(o.getXrefDataSet());
                assoc.setAssocType("weak_ortholog");
                assoc.setCreationDate(new Date());
                assoc.setMasterRgdId(o.getSrcRgdId());
                assoc.setDetailRgdId(o.getDestRgdId());
                assoc.setSrcPipeline(o.getXrefDataSrc());
                insertAssociation(assoc);

                downgradedOrthologs++;
            }

            // check if there is reverse weak ortholog (association)
            if( areGenesAssociated(o.getDestRgdId(), o.getSrcRgdId()) ) {
                reverseAssocMatch++;
                detachedAssociations += detachAssociation(forDeleteAssociations, o.getDestRgdId(), o.getSrcRgdId());
                continue; // there is reverse association!
            }

            reverseMissing++;
            Association assoc = new Association();
            assoc.setAssocSubType(o.getXrefDataSet());
            assoc.setAssocType("weak_ortholog");
            assoc.setCreationDate(new Date());
            assoc.setDetailRgdId(o.getSrcRgdId());
            assoc.setMasterRgdId(o.getDestRgdId());
            assoc.setSrcPipeline(o.getXrefDataSrc());
            insertAssociation(assoc);
        }

        deleteOrthologs(orthologsForDelete);

        String species = " ["+SpeciesType.getCommonName(speciesType1)+"==>"+SpeciesType.getCommonName(speciesType2)+"]=";
        logger.info("reverseOrthologMatch"+species+reverseOrthologMatch);
        logger.info("reverseAssocMatch   "+species+reverseAssocMatch);
        logger.info("reverseMissing      "+species+reverseMissing);
        logger.info("detachedAssociations"+species+detachedAssociations);
        logger.info("downgradedOrthologs "+species+downgradedOrthologs);
        logger.info("nondowngradedOrthos "+species+nondowngradedOrthologs);
    }

    void checkComplementAssociations(List<Association> forDeleteAssociations, int speciesTypeKey) throws Exception {

        int reverseOrthologMatch = 0;
        int reverseAssocMatch = 0;
        int reverseMissing = 0;
        int detachedAssociations = 0;

        List<Association> list = associationDAO.getAssociationsByType("weak_ortholog", SpeciesType.HUMAN, speciesTypeKey);
        for( Association a: list ) {

            // check reverse ortholog
            if( areGenesOrthologous(a.getDetailRgdId(), a.getMasterRgdId())!=0 ) {
                reverseOrthologMatch++;
                continue; // reverse ortholog does exist
            }

            // there is no reverse ortholog!
            // check if there is reverse weak ortholog (association)
            if( areGenesAssociated(a.getDetailRgdId(), a.getMasterRgdId()) ) {
                // there is also o reverse assoc in RGD;
                // if both forward and reverse association is on for-delete list, do nothing: it should be deleted
                int index1 = findAssociation(forDeleteAssociations, a.getDetailRgdId(), a.getMasterRgdId());
                int index2 = findAssociation(forDeleteAssociations, a.getMasterRgdId(), a.getDetailRgdId());
                if( index1 >=0 && index2>=0 )
                    continue; // forward and reverse assoc are both on for-delete list and they

                // either forward or reverse assoc is not on for-delete list: it should NOT be deleted
                // so we detach it from to-delete list
                reverseAssocMatch++;
                detachedAssociations += detachAssociation(forDeleteAssociations, a.getDetailRgdId(), a.getMasterRgdId());
                continue; // there is reverse association!
            }

            // there is a missing reverse association: create it
            reverseMissing++;
            Association assoc = new Association();
            assoc.setAssocSubType(a.getAssocSubType());
            assoc.setAssocType(a.getAssocType());
            assoc.setCreationDate(new Date());
            assoc.setDetailRgdId(a.getMasterRgdId());
            assoc.setMasterRgdId(a.getDetailRgdId());
            assoc.setSrcPipeline(a.getSrcPipeline());
            insertAssociation(assoc);
        }

        logger.info("reverseOrthoMatch Round2  =" +reverseOrthologMatch);
        logger.info("reverseAssocMatch Round2  =" +reverseAssocMatch);
        logger.info("reverseMissing    Round2  ="+reverseMissing);
        logger.info("detachedAssociations Round2 ="+detachedAssociations);
    }

    public void insertAssociation(Association a) throws Exception {

        validateAssociation(a);
        associationDAO.insertAssociation(a);
    }

    public void updateAssociation(Association a) throws Exception {

        validateAssociation(a);
        associationDAO.updateAssociation(a);
    }

    void validateAssociation(Association a) {
        if( !Utils.isStringEmpty(a.getAssocSubType()) ) {
            String[] sourcesOrig = a.getAssocSubType().split(", ");
            if (sourcesOrig.length > 1) {
                Set<String> sourcesUnique = new TreeSet<>(Arrays.asList(sourcesOrig));
                String newXrefDataSet = Utils.concatenate(sourcesUnique, ", ");
                if (!newXrefDataSet.equals(a.getAssocSubType()) && newXrefDataSet.length() < a.getAssocSubType().length()) {
                    a.setAssocSubType(newXrefDataSet);
                }
            }
        }
    }

    public List<Association> getAssociationsByType(String assocType, int speciesTypeKey1, int speciesTypeKey2) throws Exception {
        List<Association> assocs = associationDAO.getAssociationsByType(assocType, speciesTypeKey1, speciesTypeKey2);
        assocs.addAll(associationDAO.getAssociationsByType(assocType, speciesTypeKey2, speciesTypeKey1));
        return assocs;
    }

    public void deleteAssociationByKey(int assocKey) throws Exception {
        associationDAO.deleteAssociationByKey(assocKey);
    }

    int findAssociation(List<Association> assocs, int masterRgdId, int detailRgdId) {

        for( int i=0; i<assocs.size(); i++ ) {
            Association assoc = assocs.get(i);
            if( assoc.getMasterRgdId()==masterRgdId && assoc.getDetailRgdId()==detailRgdId ) {
                return i;
            }
        }
        return -1;
    }

    int detachAssociation(List<Association> assocs, int masterRgdId, int detailRgdId) {

        Iterator<Association> it = assocs.iterator();
        while( it.hasNext() ) {
            Association assoc = it.next();
            if( assoc.getMasterRgdId()==masterRgdId && assoc.getDetailRgdId()==detailRgdId ) {
                it.remove();
                return 1;
            }
        }
        return 0;
    }

    /**
     * get gene orthologs for given pair of species that were modified before given date
     * @param modifiedDate ortholog modification date
     * @param speciesTypeKey1 species type key1
     * @param speciesTypeKey2 species type key2
     * @return List of Ortholog objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<Ortholog> getOrthologsModifiedBefore(Date modifiedDate, int speciesTypeKey1, int speciesTypeKey2) throws Exception {
        return orthologDAO.getOrthologsModifiedBefore(modifiedDate, speciesTypeKey1, speciesTypeKey2);
    }

    /**
     * get gene orthologs for given ortholog group id
     * @param groupId group id
     * @return List of Ortholog objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<Ortholog> getOrthologsForGroupId(int groupId, String xrefDataSrc) throws Exception {
        List<Ortholog> orthologs = orthologDAO.getOrthologsForGroupId(groupId);
        // remove non-homolegene orthologs from the result list
        Iterator<Ortholog> it = orthologs.iterator();
        while( it.hasNext() ) {
            Ortholog o = it.next();
            if( !o.getXrefDataSrc().equals(xrefDataSrc) )
                it.remove();
        }
        return orthologs;
    }

    /**
     * insert a list of orthologs; Note: ortholog genetogene_key will be taken from sequence GENETOGENE_RGD_ID_RLT_SEQ
     * @param orthologs list of orthologs to insert
     * @return count of rows inserted
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int insertOrthologs(List<Ortholog> orthologs) throws Exception {
        for( Ortholog o: orthologs ) {
            // validate ortholog's xref data set
            if( !Utils.isStringEmpty(o.getXrefDataSet()) ) {
                String[] sourcesOrig = o.getXrefDataSet().split(", ");
                if (sourcesOrig.length > 1) {
                    Set<String> sourcesUnique = new TreeSet<>(Arrays.asList(sourcesOrig));
                    String newXrefDataSet = Utils.concatenate(sourcesUnique, ", ");
                    if (!newXrefDataSet.equals(o.getXrefDataSet()) && newXrefDataSet.length() < o.getXrefDataSet().length()) {
                        o.setXrefDataSet(newXrefDataSet);
                    }
                }
            }

            // set the correct ortholog type:
            //  directOrthologType for human-species relations
            //  transitiveOrthologType for non-human relations
            if( o.getSrcSpeciesTypeKey()==SpeciesType.HUMAN || o.getDestSpeciesTypeKey()==SpeciesType.HUMAN ) {
                o.setOrthologTypeKey(getDirectOrthologTypeKey());
            } else {
                o.setOrthologTypeKey(getTransitiveOrthologTypeKey());
            }

            logInsertedOrthologs.info(o.dump("|"));
        }
        return orthologDAO.insertOrthologs(orthologs);
    }

    /**
     * delete a list of orthologs; Note: ortholog key will be used
     * @param orthologs list of orthologs to  delete
     * @return count of rows deleted
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int deleteOrthologs(List<Ortholog> orthologs) throws Exception {

        // never delete custom RGD ortholog
        Iterator<Ortholog> it = orthologs.iterator();
        while( it.hasNext() ) {
            Ortholog o = it.next();
            if( Utils.stringsAreEqual(o.getXrefDataSrc(), "RGD") ) {
                it.remove();
            }
        }

        for( Ortholog o: orthologs ) {
            logDeletedOrthologs.info(o.dump("|"));
        }
        return orthologDAO.deleteOrthologs(orthologs);
    }

    /**
     * update last modified date and user (last_modified_by) for a list of orthologs
     * @param orthologs list of orthologs to be updated
     * @return count of rows updated
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int updateLastModified(List<Ortholog> orthologs, int lastModifiedBy) throws Exception {
        return orthologDAO.updateLastModified(orthologs, lastModifiedBy);
    }


    public List<Gene> getGenesByGeneId(String geneId) throws Exception {
        return getGenesByXdbId(geneId, XdbId.XDB_KEY_ENTREZGENE);
    }

    public List<Integer> getRgdIdListByEGID(String egId) throws Exception {
        JdbcTemplate jt = new JdbcTemplate(this.getDataSource());
        return jt.queryForList("select distinct x.RGD_ID from RGD_ACC_XDB x,GENES g where xdb_key=3 and ACC_ID=? "
                +"and x.RGD_ID=g.RGD_ID and g.GENE_TYPE_LC NOT IN('allele','splice')",
                new Object[]{egId}, new int[]{Types.VARCHAR}, Integer.class);
    }

    public List<Gene> getGenesByXdbId(String geneId, int xdbKey) throws Exception {
        String key = xdbKey + geneId;
        List<Gene> genes = _genesByXdbIdCache.get(key);
        if( genes==null ) {
            genes = xdbIdDAO.getActiveGenesByXdbId(xdbKey, geneId);
            _genesByXdbIdCache.put(key, genes);
        }
        return genes;
    }
    private Map<String, List<Gene>> _genesByXdbIdCache = new HashMap<>();


    /**
     * get status of object identified by rgd id
     * @param rgdId rgd id of object
     * @return object status: one of ('ACTIVE','WITHDRAWN','RETIRED')
     * @throws Exception thrown if rgd id is invalid
     */
    public String getStatus(int rgdId) throws Exception {
        RgdId rgdIdObj = getRgdId(rgdId);
        return rgdIdObj==null ? null : rgdIdObj.getObjectStatus();
    }

    synchronized RgdId getRgdId(int rgdId) throws Exception {

        RgdId rgdIdObj = _rgdIdCache.get(rgdId);
        if( rgdIdObj==null ) {
            rgdIdObj = rgdIdDAO.getRgdId2(rgdId);
            _rgdIdCache.put(rgdId, rgdIdObj);
        }
        return rgdIdObj;
    }
    static Map<Integer, RgdId> _rgdIdCache = new HashMap<>(10003);

    Gene getGeneByRgdId(int rgdId) throws Exception {
        return geneDAO.getGene(rgdId);
    }

    synchronized public String getGeneSymbol(int geneRgdId) {

        String geneSymbol = _geneSymbolCache.get(geneRgdId);
        if( geneSymbol==null ) {
            Gene gene = null;
            try {gene = geneDAO.getGene(geneRgdId); }
            catch(Exception ignore) {}
            if( gene!=null ) {
                geneSymbol = gene.getSymbol();
                _geneSymbolCache.put(geneRgdId, geneSymbol);
            }
        }
        return geneSymbol;
    }
    static Map<Integer, String> _geneSymbolCache = new HashMap<>(10003);

    public Gene getGeneBySymbol(String geneSymbol, int speciesTypeKey, Logger log) throws Exception {
        List<Gene> genes = geneDAO.getAllGenesBySymbol(geneSymbol, speciesTypeKey);
        if( genes.size()>1 ) {
            // remove inactive genes
            Iterator<Gene> it = genes.iterator();
            while( it.hasNext() ) {
                Gene g = it.next();
                RgdId id = getRgdId(g.getRgdId());
                if( !id.getObjectStatus().equals("ACTIVE") ) {
                    it.remove();
                }
            }
        }
        if( genes.size()>1 ) {
            log.warn("multiple genes for symbol "+geneSymbol+", species "+speciesTypeKey);
            return null;
        }
        return genes.isEmpty() ? null : genes.get(0);
    }

    public Gene insertAgrGene(int speciesTypeKey, String geneSymbol, String agrGeneId) throws Exception {

        // extra security: cannot insert a HUMAN, MOUSE, RAT gene
        if( speciesTypeKey==SpeciesType.HUMAN || speciesTypeKey==SpeciesType.MOUSE || speciesTypeKey==SpeciesType.RAT ) {
            return null;
        }

        RgdId id = rgdIdDAO.createRgdId(RgdId.OBJECT_KEY_GENES, "ACTIVE", "created by AGR Ortholog Loader", speciesTypeKey);

        Gene gene = new Gene();
        gene.setRgdId(id.getRgdId());
        gene.setSymbol(geneSymbol);
        gene.setSpeciesTypeKey(speciesTypeKey);
        geneDAO.insertGene(gene);

        insertAgrGeneXdbId(id.getRgdId(), agrGeneId);

        Logger log = Logger.getLogger("insertedAgrGenes");
        log.debug("RGD:"+id.getRgdId()+" SYMBOL:"+geneSymbol+" AGR_CURIE:"+agrGeneId);

        return gene;
    }

    void insertAgrGeneXdbId(int geneRgdId, String accId) throws Exception {
        XdbId xdbId = new XdbId();
        xdbId.setRgdId(geneRgdId);
        xdbId.setAccId(accId);
        xdbId.setXdbKey(63);
        xdbId.setSrcPipeline("AgrOrtholog");
        xdbIdDAO.insertXdb(xdbId);
    }

    public List<String> getCrossLinkedOrthologs(int speciesTypeKey) throws Exception {
        JdbcTemplate jt=new JdbcTemplate(this.getDataSource());
        return jt.query("SELECT g.rgd_id,g.gene_symbol FROM genes g WHERE g.rgd_id IN(\n" +
                "SELECT dest_rgd_id FROM GENETOGENE_RGD_ID_RLT,rgd_ids r1,rgd_ids r2\n" +
                "WHERE dest_rgd_id=r1.rgd_id AND r1.object_status='ACTIVE' AND r1.species_type_key=?\n" +
                "  AND src_rgd_id=r2.rgd_id AND r2.object_status='ACTIVE' AND r2.species_type_key IN(1,2,3)\n" +
                "GROUP BY dest_rgd_id HAVING COUNT(*)>2\n" +
                "UNION\n" +
                "SELECT src_rgd_id FROM GENETOGENE_RGD_ID_RLT,rgd_ids r1,rgd_ids r2\n" +
                "WHERE src_rgd_id=r1.rgd_id AND r1.object_status='ACTIVE' AND r1.species_type_key=?\n" +
                "  AND dest_rgd_id=r2.rgd_id AND r2.object_status='ACTIVE' AND r2.species_type_key IN(1,2,3)\n" +
                "GROUP BY src_rgd_id HAVING COUNT(*)>2\n" +
                ") ORDER BY LOWER(gene_symbol)",
                new Object[]{speciesTypeKey, speciesTypeKey},
                new RowMapper(){
                    public Object mapRow(ResultSet resultSet, int i) throws SQLException {
                        return resultSet.getInt(1)+"\t"+resultSet.getString(2);
                    }
                });
    }

    /**
     * get all active genes without any orthologs for given species
     * Note: splices and alleles are excluded from the list
     * @param speciesTypeKey species type key
     * @return List of Gene objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<Gene> getGenesWithoutOrthologs(int speciesTypeKey) throws Exception {
        return orthologDAO.getGenesWithoutOrthologs(speciesTypeKey);
    }

    public int getReplacedGeneByRgdId(int oldRgdId) throws Exception {

        return rgdIdDAO.getActiveRgdIdFromHistory(oldRgdId);
    }

    public int fixXrefDataSetInOrthologs() throws Exception {

        int rowsProcessed = 0;
        int rowsUpdated = 0;

        List<Ortholog> orthologs = orthologDAO.getOrthologsModifiedBefore(new Date());
        for( Ortholog o: orthologs ) {
            rowsProcessed++;

            if( Utils.isStringEmpty(o.getXrefDataSet()) ) {
                continue;
            }

            // compute non-duplicate list of sources from the current XREF_DATA_SET
            String[] sourcesOrig = o.getXrefDataSet().split(", ");
            if( sourcesOrig.length==1 ) {
                continue;
            }
            Set<String> sourcesUnique = new TreeSet<>(Arrays.asList(sourcesOrig));
            String newXrefDataSet = Utils.concatenate(sourcesUnique, ", ");
            if( !newXrefDataSet.equals(o.getXrefDataSet()) && newXrefDataSet.length()<o.getXrefDataSet().length() ) {
                String sql = "UPDATE genetogene_rgd_id_rlt SET xref_data_set=? WHERE genetogene_key=?";
                orthologDAO.update(sql, newXrefDataSet, o.getKey());
                rowsUpdated++;
            }
        }

        logger.info("orthologs processed: "+rowsProcessed);
        return rowsUpdated;
    }

    public int fixXrefDataSetInAssociations() throws Exception {

        int rowsProcessed = 0;
        int rowsUpdated = 0;

        List<Association> assocs = associationDAO.getAssociationsByType("weak_ortholog");
        for( Association a: assocs ) {
            rowsProcessed++;

            if( Utils.isStringEmpty(a.getAssocSubType()) ) {
                continue;
            }

            // compute non-duplicate list of sources from the current XREF_DATA_SET
            String[] sourcesOrig = a.getAssocSubType().split(", ");
            if( sourcesOrig.length==1 ) {
                continue;
            }
            Set<String> sourcesUnique = new TreeSet<>(Arrays.asList(sourcesOrig));
            String newXrefDataSet = Utils.concatenate(sourcesUnique, ", ");
            if( !newXrefDataSet.equals(a.getAssocSubType()) && newXrefDataSet.length()<a.getAssocSubType().length() ) {
                a.setAssocSubType(newXrefDataSet);
                associationDAO.updateAssociation(a);
                rowsUpdated++;
            }
        }

        logger.info("associations processed: "+rowsProcessed);
        return rowsUpdated;
    }

    /**
     * get ortholog count for given pair of species
     * @param speciesTypeKey1 species type key for first species
     * @param speciesTypeKey2 species type key for second species
     * @return count of orthologs
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int getOrthologCount(int speciesTypeKey1, int speciesTypeKey2) throws Exception {

        return orthologDAO.getOrthologCount(speciesTypeKey1, speciesTypeKey2);
    }

    /**
     * return count of associations for given assoc type and pair of species
     * @param assocType assoc type
     * @param masterSpeciesTypeKey species type key for master rgd id
     * @param detailSpeciesTypeKey species type key for detail rgd id
     * @return count of associations for given assoc type and pair of species
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int getAssociationCountByType(String assocType, int masterSpeciesTypeKey, int detailSpeciesTypeKey) throws Exception {

        return associationDAO.getAssociationCountByType(assocType, masterSpeciesTypeKey, detailSpeciesTypeKey);
    }

    public DataSource getDataSource() throws Exception {
        return orthologDAO.getDataSource();
    }

    public void setDirectOrthologTypeKey(int directOrthologTypeKey) {
        this.directOrthologTypeKey = directOrthologTypeKey;
    }

    public int getDirectOrthologTypeKey() {
        return directOrthologTypeKey;
    }

    public void setTransitiveOrthologTypeKey(int transitiveOrthologTypeKey) {
        this.transitiveOrthologTypeKey = transitiveOrthologTypeKey;
    }

    public int getTransitiveOrthologTypeKey() {
        return transitiveOrthologTypeKey;
    }
}
