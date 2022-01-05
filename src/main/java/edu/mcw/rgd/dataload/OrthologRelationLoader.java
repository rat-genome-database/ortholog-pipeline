package edu.mcw.rgd.dataload;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import edu.mcw.rgd.datamodel.Association;
import edu.mcw.rgd.datamodel.Ortholog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/*
 * Created on Mar 23, 2007
 * @author dli
 * create the homolog group id, match the gene in the homolog to RGD entries
 */
public class OrthologRelationLoader {
    OrthologRelationParser orthologRelationParser;
    OrthologRelationDao dao;
    List<OrthologRelation> orthologRelations;
    final String unmatched="unmatched";
    final String multiple="multiple";
    final String withdrawn="withdrawn";
    CounterPool counters;

    protected final Logger loggerMatched = LogManager.getLogger("matched");
    protected final Logger loggerUnmatch = LogManager.getLogger("unmatched");
    protected final Logger multipleMatch = LogManager.getLogger("multipleMatch");
    protected final Logger loggerWithdrawn = LogManager.getLogger("withdrawn");
    protected final Logger process = LogManager.getLogger("process");
    protected final Logger assoc = LogManager.getLogger("assoc");
    private String defaultDataSetName;
    private int createdBy; // curation user id for the pipeline

    public void run(List<OrthologRelation> orthologRelations, int speciesTypeKey) throws Exception {

        counters = new CounterPool();

        setOrthologRelations(orthologRelations);

        matchRgdId();
        loadData(speciesTypeKey);
    }

    // updated logic:
    // OrthologGroup is grouped by human rgd id; it contains one ortholog per species
    // the other orthologs are weak ortholog
    // ortholog entries within a group come from RGD manual orthologs and incoming data
    // for every species one best consensus entry is chosen as ortholog, the other go as weak orthologs
    public void loadData(int speciesTypeKey) throws Exception {

        // drop unmapped orthologs from initial ortholog relations collection
        dropUnmappedRelations();

        // build ortholog groups -- group id is a human gene rgd id
        // only best-fit relationship becomes an ortholog, a strong ortholog
        // the rest becomes other ortholog, or weak ortholog
        // also note that RGD manual ortholog suppresses any other orthologs

        Collection<OrthologGroup> groups = buildGroups(orthologRelations);

        qcGroups(groups, speciesTypeKey);

        List<Ortholog> weakOrthologs = loadGroups(groups, speciesTypeKey);

        // extract non-orthologous relations and turn them into associations
        List<Association> associations = buildAssociationsFromRelations(groups);
        addWeakOrthologs(associations, weakOrthologs);
        dropStrongOrthologsFromAssociations(associations);
        processOrthologAssociations(associations, speciesTypeKey);
    }

    // drop relations that could not be mapped to RGD
    void dropUnmappedRelations() throws Exception {

        int countOriginal = orthologRelations.size();
        Iterator<OrthologRelation> it = orthologRelations.iterator();
        while( it.hasNext() ) {
            OrthologRelation orRel = it.next();

            if( orRel.getDestRgdId()<=0 || orRel.getSrcRgdId()<=0 ) {
                it.remove();
            }
        }
        int unmappedRelations = countOriginal - orthologRelations.size();
        process.info("Dropped unmapped relations: "+unmappedRelations);
    }

    void addWeakOrthologs(List<Association> associations, List<Ortholog> weakOrthologs) {

        process.info("  building associations from weak orthologs: "+weakOrthologs.size());

        for( Ortholog o: weakOrthologs ) {
            Association assoc = new Association();
            assoc.setAssocSubType(o.getXrefDataSet());
            assoc.setAssocType("weak_ortholog");
            assoc.setCreationDate(o.getCreatedDate());
            assoc.setDetailRgdId(o.getDestRgdId());
            assoc.setMasterRgdId(o.getSrcRgdId());
            assoc.setSrcPipeline(o.getXrefDataSrc());
            associations.add(assoc);
        }

        process.info("  INCOMING ASSOCIATIONS + WEAK ORTHOLOGS: " + associations.size());
    }

    List<Association> buildAssociationsFromRelations(Collection<OrthologGroup> groups) throws Exception {

        List<Association> associations = new ArrayList<>(orthologRelations.size());

        for( OrthologGroup group: groups ) {
            for( OrthologRelation orRel: group.relations ) {
                Association assoc = new Association();
                assoc.setAssocType("weak_ortholog");
                assoc.setAssocSubType(orRel.getDataSetName());
                assoc.setSrcPipeline(orRel.getDataSource());
                assoc.setMasterRgdId(orRel.getSrcRgdId());
                assoc.setDetailRgdId(orRel.getDestRgdId());
                assoc.setCreationDate(new Date());
                associations.add(assoc);
            }
        }

        process.info("  INCOMING ASSOCIATIONS: " + associations.size());

        return associations;
    }

    void dropStrongOrthologsFromAssociations(List<Association> associations) throws Exception {

        Iterator<Association> it = associations.iterator();
        while( it.hasNext() ) {
            Association assoc = it.next();
            int orthologyKey = dao.areGenesOrthologous(assoc.getMasterRgdId(), assoc.getDetailRgdId());
            if( orthologyKey!=0 ) {
                it.remove();
            }
        }
        process.info("  INCOMING ASSOCIATIONS + WEAK ORTHOS - STRONG ORTHOS: "+associations.size());
    }

    void processOrthologAssociations(List<Association> associations, int speciesTypeKey) throws Exception {

        List<Association> assocInRgd = dao.getAssociationsByType("weak_ortholog", speciesTypeKey, SpeciesType.HUMAN);

        OrthoAssociationSyncer syncer = new OrthoAssociationSyncer();
        syncer.addIncomingObjects(associations);
        syncer.setObjectsInRgd(assocInRgd);

        // run qc - delete stale objects afterwards (in-rgd associations that have not been processed by qc process)
        boolean deleteStaleObjects = true;
        syncer.qc("assoc", deleteStaleObjects);

        // insert map data
        List<Association> objs = syncer.getObjectsForInsert();
        if( !objs.isEmpty() ) {
            for( Association assoc: objs ) {
                this.assoc.info("INSERT|"+assoc.dump("|"));
                dao.insertAssociation(assoc);
            }
            process.info("WEAK_ORTHOLOGS_INSERTED: " + objs.size());
        }
        // update map data
        objs = syncer.getObjectsForUpdate();
        if( !objs.isEmpty() ) {
            for( Association assoc: objs ) {
                this.assoc.info("UPDATE|"+assoc.dump("|"));
                dao.updateAssociation(assoc);
            }
            process.info("WEAK_ORTHOLOGS_UPDATED: "+objs.size());
        }

        objs = syncer.getObjectsForDelete();

        // we run detailed QC checks for all orthologs; if there is any reverse association matching association to delete
        // then we found a NO-OP: deletion followed by insertion; we remove these associations from to-delete list
        dao.checkComplementOrthologs(SpeciesType.HUMAN, speciesTypeKey, objs);
        dao.checkComplementOrthologs(speciesTypeKey, SpeciesType.HUMAN, objs);

        // we run detailed QC checks for all association; if there is any reverse association matching association to delete
        // then we found a NO-OP: deletion followed by insertion; we remove these associations from to-delete list
        dao.checkComplementAssociations(objs, speciesTypeKey);

        // delete map data
        if( !objs.isEmpty() ) {

            for( Association assoc: objs ) {
                this.assoc.info("DELETE|"+assoc.dump("|"));
                dao.deleteAssociationByKey(assoc.getAssocKey());
            }
            process.info("WEAK_ORTHOLOGS_DELETED: " + objs.size());
        }

        process.info("WEAK_ORTHOLOGS_MATCHED: " + syncer.getMatchingObjects().size());
    }

    /*
     * map the entrezgene id in the ortholog relationship data to RGD ID
     */
    public void matchRgdId() throws Exception {

        process.info("Match relations against RGD:");

        // must be synchronized because multiple threads could modify the map at the same time
        final Map<String, String> egIdRgdId = new ConcurrentHashMap<>();

        List<OrthologRelation> incomingRecords = new ArrayList<>(orthologRelations.size());
        for (OrthologRelation orRel:orthologRelations) {
            incomingRecords.add(orRel);
        }

        incomingRecords.parallelStream().forEach( orRel -> {

            String idSrc=orRel.getSrcOtherId(); //the source eg id
            String idDest=orRel.getDestOtherId(); // the destination eg id
            String srcRgdId; // the source rgd id
            String destRgdId; // the source eg id

            process.debug("QC EGSRC="+idSrc+", EGDST="+idDest);

            // processing source first, if the source couldn't be matched, don't process destination any more
            // first check if the query is already made and stored in egIdRgdId map
            if (egIdRgdId.get(idSrc)==null) {
                // if the query is not made yet, make the query and store it in the map
                String rgdId=getRgdIdByEgId(idSrc);
                egIdRgdId.put(idSrc,rgdId);
                srcRgdId=rgdId;
            } else {
                // the query is already made, take it from the map
                srcRgdId = egIdRgdId.get(idSrc);
            }
            if (srcRgdId.contains(multiple)) {
                multipleMatch.error("EGID:"+idSrc+"\t"+srcRgdId+"\t"+orRel.getSrcSpeciesTypeKey());
                counters.increment("multipleMatchC");
                return;
            }
            else if (srcRgdId.contains(unmatched)) {
                loggerUnmatch.error("EGID:"+ idSrc +"\t" +orRel.getSrcSpeciesTypeKey());
                counters.increment("unMatchedC");
                return;
            }
            else if (srcRgdId.contains(withdrawn)) {
                loggerWithdrawn.error("EGID:"+ idSrc +"\t"+srcRgdId+"\t" +orRel.getSrcSpeciesTypeKey());
                counters.increment("withdrawnC");
                return;
            }
            //processing destination
            // first check if the query is already made and stored in egIdRgdId map
            if (egIdRgdId.get(idDest)==null) {
                // if the query is not made yet, make the query and store it in the map
                String rgdId=getRgdIdByEgId(idDest);
                egIdRgdId.put(idDest,rgdId);
                destRgdId=rgdId;
            } else {
                // the query is already made, take it from the map
                destRgdId = egIdRgdId.get(idDest);
            }

            if (destRgdId.contains(multiple)) {
                multipleMatch.error("EGID:"+idDest+"\t"+destRgdId + "\t"+ orRel.getDestSpeciesTypeKey());
                counters.increment("multipleMatchC");
                return;
            } else if (destRgdId.contains(unmatched)) {
                loggerUnmatch.error("EGID:"+ idDest + "\t"+ orRel.getDestSpeciesTypeKey());
                counters.increment("unMatchedC");
                return;
            } else if (destRgdId.contains(withdrawn)) {
                loggerWithdrawn.error("EGID:"+ idDest +"\t"+destRgdId+"\t" +orRel.getDestSpeciesTypeKey());
                counters.increment("withdrawnC");
                return;
            }
            orRel.setSrcRgdId(Integer.parseInt(srcRgdId));
            orRel.setDestRgdId(Integer.parseInt(destRgdId));
            loggerMatched.info("src EGID:"+idSrc+", src RGD:"+srcRgdId + ", dest EGID:"+idDest+", dest RGD:"+destRgdId +", src species:"+orRel.getSrcSpeciesTypeKey()+", dest species:"+orRel.getDestSpeciesTypeKey());
            counters.increment("matchedC");
        });

        process.info("  matched relations:       " + getMatchedC());
        process.info("  unmatched relations:     " + getUnMatchedC());
        process.info("  matched withdrawn genes: " + getWithdrawnC());
        process.info("  multiple match relations:" + getMultipleMatchC());
    }

    private String getRgdIdByEgId(String egId) {
        try {
            return getRgdIdByEgId2(egId);
        } catch( Exception e ) {
            throw new RuntimeException(e);
        }
    }

    /* In: an entrez gene id
     * Returns: 
     *   unmatched ----- the entrezgene id doesn't matches any genes in RGD (including active and non active genes)
     *   multiple ---- the entrezgene id matches multiple genes( not including splice and allele genes) or matche multiple genes
     *                 that are all splice or allele genes
     *   withdrawn ---- the entrezgene id matches non active gene which is not replaced by any active gene
     *   rgd id ----- the entrezgene id matches only one active gene or one nonactive gene which has been replaced by an active gene
     */
    private String getRgdIdByEgId2(String egId) throws Exception {
        
        String rRgdId=""; // the returned rgd id
        String rgdIdInfo=""; // the associated rgd ids if there is multiple match
        
        // the rgdIds queried by the entrezgene id; splice and allele are excluded
        List<Integer> rgdIds= dao.getRgdIdListByEGID(egId);
        //logger.debug("rgd ids: "+ rgdIds);
        List<Integer> nonActiveRgdIds= new ArrayList<>();
        
        if (rgdIds==null || rgdIds.isEmpty() )
            return unmatched;

        // the result has rgd id
        int activeC=0;     // count active genes
        for (int rgdId: rgdIds ) {
            // get status for each gene
            String status= dao.getStatus(rgdId);
            if (!status.equals("ACTIVE")) {
                nonActiveRgdIds.add(rgdId);
                continue;
            }
            rRgdId=rgdId+"";
            activeC++;
            rgdIdInfo=rgdIdInfo+","+rRgdId;
        }

        if ( activeC >1) {
            // the gene is associated with multiple active genes
            return multiple+":"+rgdIdInfo;
        } else if (activeC ==1) {
            // the gene matches only one active gene
            return rRgdId;
        } else {
            // activeC=0, the gene is associated with only nonactive genes
            int replacedRgdId=0;
            int replacedActiveC=0;
            int nonActiveRgdId=0;
            String reRgdIdInfo="";

            // find the replaced genes
            List<Integer> replacedRgdIds= new ArrayList<>();
            for (int i=0; i<nonActiveRgdIds.size(); i++) {
                nonActiveRgdId=nonActiveRgdIds.get(i);
                replacedRgdId= dao.getReplacedGeneByRgdId(nonActiveRgdId);
                if (replacedRgdId >0) {
                    replacedRgdIds.add(replacedRgdId); // the replaced gene is active
                    replacedActiveC++;
                    reRgdIdInfo=reRgdIdInfo+","+replacedRgdId;
                }
            }

            if (replacedRgdIds.size() >1) {
                    // multiple replaced genes
                    return multiple+" replaced:"+ reRgdIdInfo;
            } else if (replacedRgdIds.size() ==1) {
                // matches one replaced gene
                    return replacedRgdId +"";
            } else {
                // filteredReplacedRgdIds.size()=0, couldn't find replaced gene
                return withdrawn+":"+nonActiveRgdId;
            }
        }
    }


    Collection<OrthologGroup> buildGroups(Collection<OrthologRelation> relations) {

        // combine files in ortholog groups: group is identified by human rgd id
        int duplicateRelations = 0;
        Map<Integer, OrthologGroup> groups = new HashMap<>();
        for( OrthologRelation rel: relations ) {

            int groupId = rel.getSrcRgdId();
            if( rel.getSrcSpeciesTypeKey()!=1 )
                throw new RuntimeException("unexpected species type key");

            OrthologGroup group = groups.get(groupId);
            if( group==null ) {
                group = new OrthologGroup();
                groups.put(groupId, group);
            }
            if( !group.add(rel) ) {
                duplicateRelations++;
            }
        }

        process.info("ortholog groups created: " + groups.size()+", duplicate relations merged: "+duplicateRelations);
        return groups.values();
    }

    // 1. separate group relations by species
    // 2. for every species subgroup, choose an ortholog
    // 3. remaining entries for species subgroup will become associations
    void qcGroups( Collection<OrthologGroup> groups, int speciesTypeKey ) throws Exception {

        int[] methodCounters = new int[4];
        for( OrthologGroup group: groups ) {

            // build species-human complementary relations
            group.buildComplementaryRelations();

            // extract relations by species pairs: human-rat, etc
            // species pair is 'src_species_type_key|dest_species_type_key'
            Map<String, List<OrthologRelation>> speciesMap = extractRelationsBySpeciesPairs(group);

            // pick an ortholog for every relation group in a species pair
            for( Map.Entry<String, List<OrthologRelation>> entry: speciesMap.entrySet() ) {
                List<OrthologRelation> relations = entry.getValue();
                Ortholog ortholog = generateOrtholog(relations, methodCounters);
                if( ortholog!=null )
                    group.incomingList.add(ortholog);
            }
        }

        process.info("ORTHOLOG BEST FIT STATS:");
        process.info("  bestFitOneRel " + methodCounters[0]);
        process.info("  bestFitLongestEvidence " + methodCounters[1]);
        process.info("  bestFitSymbolMatch " + methodCounters[2]);
        process.info("  bestFitShortestSymbol " + methodCounters[3]);
    }

    // extract relations by species pairs
    Map<String, List<OrthologRelation>> extractRelationsBySpeciesPairs(OrthologGroup group) {

        Map<String, List<OrthologRelation>> speciesMap = new HashMap<>();
        for( OrthologRelation rel: group.relations ) {
            String key = rel.getSrcSpeciesTypeKey()+"|"+rel.getDestSpeciesTypeKey();
            List<OrthologRelation> relations = speciesMap.get(key);
            if( relations==null ) {
                relations = new ArrayList<>();
                speciesMap.put(key, relations);
            }
            relations.add(rel);
        }
        return speciesMap;
    }

    Ortholog generateOrtholog(List<OrthologRelation> relations, int[] methodCounters) throws Exception {

        int srcRgdId = relations.get(0).getSrcRgdId();
        int destSpeciesTypeKey = relations.get(0).getDestSpeciesTypeKey();

        List<Ortholog> manualOrthologs = dao.getManualOrthologs(srcRgdId, destSpeciesTypeKey);
        if( !manualOrthologs.isEmpty() ) {
            // manual orthologs always take precedence over incoming relations
            if( manualOrthologs.size()>1 ) {
                process.warn("CONFLICT: multiple manual orthologs for src_rgd_id="+srcRgdId+" and species type key "+destSpeciesTypeKey);
                return null;
            }
            return manualOrthologs.get(0);
        }

        // no manual orthologs present -- pick a best fit relation and remove that relation from relations list
        OrthologRelation bestFitRel = pickBestFitRelation(relations, methodCounters, dao);

        // turn the relation into an ortholog
        Ortholog o = new Ortholog();
        o.setXrefDataSrc(bestFitRel.getDataSource());
        o.setCreatedBy(70);
        o.setLastModifiedBy(70);
        o.setCreatedDate(new Date());
        o.setLastModifiedDate(new Date());

        o.setSrcRgdId(bestFitRel.getSrcRgdId());
        o.setDestRgdId(bestFitRel.getDestRgdId());
        o.setXrefDataSet(bestFitRel.getDataSetName());
        o.setSrcSpeciesTypeKey(bestFitRel.getSrcSpeciesTypeKey());
        o.setDestSpeciesTypeKey(bestFitRel.getDestSpeciesTypeKey());
        return o;
    }

    // methodCounters - array of how many times given best-fit method was used
    public static OrthologRelation pickBestFitRelation(List<OrthologRelation> relations, int[] methodCounters, final OrthologRelationDao dao) throws Exception {

        int sourceRgdId = relations.get(0).getSrcRgdId();

        // if there is only one relation, pick it!
        if( relations.size()==1 ) {
            methodCounters[0]++;
            return relations.remove(0);
        }

        // pick relation with longest evidence
        Collections.sort(relations, new Comparator<OrthologRelation>() {
            public int compare(OrthologRelation o1, OrthologRelation o2) {
                return getEvidenceCount(o2)-getEvidenceCount(o1);
            }
        });
        if( getEvidenceCount(relations.get(0)) > getEvidenceCount(relations.get(1)) ) {
            methodCounters[1]++;
            return relations.remove(0);
        }

        // pick relation with matching symbol
        String sourceGeneSymbol = dao.getGeneSymbol(sourceRgdId);
        for( int i=0; i<relations.size(); i++ ) {
            OrthologRelation rel = relations.get(i);
            if( Utils.stringsCompareToIgnoreCase(sourceGeneSymbol, dao.getGeneSymbol(rel.getDestRgdId()))==0 ) {
                methodCounters[2]++;
                return relations.remove(i);
            }
        }

        // sort relations by case insensitive gene symbols
        Collections.sort(relations, new Comparator<OrthologRelation>() {
            public int compare(OrthologRelation o1, OrthologRelation o2) {
                return Utils.stringsCompareToIgnoreCase(
                        dao.getGeneSymbol(o1.getDestRgdId()),
                        dao.getGeneSymbol(o2.getDestRgdId()));
            }
        });
        methodCounters[3]++;
        return relations.remove(0);
    }

    static int getEvidenceCount(OrthologRelation rel) {

        int count = 1;
        String evidences = rel.getDataSetName();
        for( int pos=0; (pos=evidences.indexOf(',', pos))>0; pos++ ) {
            count++;
        }
        return count;
    }


    List<Ortholog> loadGroups(Collection<OrthologGroup> groups, int speciesTypeKey) throws Exception {

        Date modifiedDate = new Date();

        int orthologsMatched = 0;
        int orthologsInserted = 0;
        int orthologsDeleted = 0;

        List<Ortholog> matchList = new ArrayList<>();
        List<Ortholog> insertList = new ArrayList<>();
        List<Ortholog> deleteList = new ArrayList<>();
        List<Ortholog> weakOrthologs = new ArrayList<>();

        dao.logDeletedOrthologs.info("LOAD GROUPS STARTING");

        for( OrthologGroup group: groups ) {

            for( Ortholog o: group.incomingList ) {

                int genetogeneKey = dao.getKeyForMatchingOrtholog(o, deleteList);
                if( genetogeneKey>0 ) {
                    // there is a matching ortholog
                    o.setKey(genetogeneKey);
                    matchList.add(o);
                }
                else if( genetogeneKey==0 ) {
                    // there is no matching ortholog and no ortholog for srcRgdId and destSpeciesTypeKey
                    // it is safe to add the ortholog
                    insertList.add(o);
                }
                else if( genetogeneKey<0 ) {
                    // there is no matching ortholog but there is annother ortholog for same srcRgdId and destSpeciesTypeKey
                    // ortholog in question must be downgraded to weak ortholog
                    weakOrthologs.add(o);
                }
            }

            if( !matchList.isEmpty() ) {
                dao.updateLastModified(matchList, 70);
                orthologsMatched += matchList.size();
                matchList.clear();
            }

            if( !insertList.isEmpty() ) {
                dao.insertOrthologs(insertList);
                orthologsInserted += insertList.size();
                insertList.clear();
            }

            if( !deleteList.isEmpty() ) {
                dao.deleteOrthologs(deleteList);
                orthologsDeleted += deleteList.size();
                deleteList.clear();
            }
        }
        dao.logDeletedOrthologs.info("LOAD GROUPS ENDED");

        process.info("ORTHOLOGS MATCHED  " + orthologsMatched);
        process.info("ORTHOLOGS INSERTED " + orthologsInserted);
        process.info("ORTHOLOGS DELETED  " + orthologsDeleted);

        handleStaleOrthologs(modifiedDate, speciesTypeKey);

        return weakOrthologs;
    }

    void handleStaleOrthologs(Date modifiedDate, int speciesTypeKey) throws Exception {

        List<Ortholog> deleteList = dao.getOrthologsModifiedBefore(modifiedDate, speciesTypeKey, SpeciesType.HUMAN);
        if( !deleteList.isEmpty() ) {
            int staleOrthologsDeleted = 0;
            int staleOrthologsSkipped = 0;
            for( Ortholog o: deleteList ) {
                if( dao.deleteStaleOrtholog(o) )
                    staleOrthologsDeleted++;
                else
                    staleOrthologsSkipped++;
            }
            process.info("STALE ORTHOLOGS DELETED " + staleOrthologsDeleted);
            process.info("STALE ORTHOLOGS SKIPPED " + staleOrthologsSkipped);
        }
    }

    public List<OrthologRelation> getOrthologRelations() {
        return orthologRelations;
    }
    public void setOrthologRelations(List<OrthologRelation> orthologRelations) {
        this.orthologRelations = orthologRelations;
    }
    public OrthologRelationParser getOrthologRelationParser() {
        return orthologRelationParser;
    }
    public void setOrthologRelationParser(OrthologRelationParser orthologRelationParser) {
        this.orthologRelationParser = orthologRelationParser;
    }
    
    public OrthologRelationDao getOrthologRelationDao() {
        return dao;
    }
    public void setOrthologRelationDao(OrthologRelationDao orthologRelationDao) {
        this.dao = orthologRelationDao;
    }
    
    
    public int getMatchedC() {
        return counters.get("matchedC");
    }

    public int getMultipleMatchC() {
        return counters.get("multipleMatchC");
    }

    public int getUnMatchedC() {
        return counters.get("unMatchedC");
    }

    public int getWithdrawnC() {
        return counters.get("withdrawnC");
    }
    
    public int getRowsDeleted() {
        return counters.get("rowsDeleted");
    }

    public int getRowsInserted() {
        return counters.get("rowsInserted");
    }

    public void setDefaultDataSetName(String defaultDataSetName) {
        this.defaultDataSetName = defaultDataSetName;
    }

    public String getDefaultDataSetName() {
        return defaultDataSetName;
    }

    public int getCreatedBy() {
        return createdBy;
    }

    public void setCreatedBy(int createdBy) {
        this.createdBy = createdBy;
    }

}
