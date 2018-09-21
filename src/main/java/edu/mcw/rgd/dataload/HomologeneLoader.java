package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Ortholog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.sync.AssociationSyncer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 5/13/13
 * Time: 4:29 PM
 * load ortholog entries from homologene
 * <p>
 *     Note: not used on production anymore
 * </p>
 */
@Deprecated
public class HomologeneLoader {

    OrthologRelationDao dao = new OrthologRelationDao();

    int badSpecies; // how many homolegen entries skipped because of unsupported species
    int noRgdId; // count of homologene entries no matching rgd id
    int multiRgdId; // count of homologene entries with multiple matching rgd id

    public static void main(String[] args) throws Exception {

        HomologeneLoader loader = new HomologeneLoader();
        loader.run();
    }

    void run() throws Exception {

        // download source file
        String localFile = downloadSourceFile();

        // parse file
        List<HOrthologEntry> entries = parseFile(localFile);

        // combine files in ortholog groups
        Collection<HOrthologGroup> groups = buildGroups(entries);

        qcGroups(groups);

        loadGroups(groups);

        System.out.println("OK");
    }

    String downloadSourceFile() throws Exception {

        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data");
        downloader.setAppendDateStamp(true);
        downloader.setLocalFile("data/homologene.data");

        System.out.println("Downloading file " + downloader.getExternalFile());
        String localFile = downloader.download();
        System.out.println("File has been downloaded to "+localFile);

        return localFile;
    }

    List<HOrthologEntry> parseFile(String fileName) throws Exception {

        List<HOrthologEntry> entries = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line;
        while ((line=br.readLine())!=null) {

            String[] cols = line.split("[\t]", -1);

            HOrthologEntry entry = new HOrthologEntry();
            entry.setHomologeneId(Integer.parseInt(cols[0]));
            entry.setTaxonomicId(Integer.parseInt(cols[1]));
            entry.setGeneId(cols[2]);
            entry.setDataSet("HOMOLOGENE");

            if( entry.getSpeciesTypeKey()== SpeciesType.ALL )
                badSpecies++;
            else
                entries.add(entry);
        }
        br.close();
        System.out.println(entries.size() + " good entries were read from homologene.data file");
        System.out.println(badSpecies+" bad species entries were read from homologene.data file");

        return entries;
    }

    Collection<HOrthologGroup> buildGroups(Collection<HOrthologEntry> entries) {

        // combine files in ortholog groups
        Map<Integer, HOrthologGroup> groups = new HashMap<>();
        for( HOrthologEntry entry: entries ) {

            int groupId = entry.getHomologeneId();
            HOrthologGroup group = groups.get(groupId);
            if( group==null ) {
                group = new HOrthologGroup();
                group.setHomolegeneId(groupId);
                groups.put(groupId, group);
            }
            group.getEntries().add(entry);
        }
        return groups.values();
    }

    void qcGroups(Collection<HOrthologGroup> groups) throws Exception {

        for( HOrthologGroup group: groups ) {

            // match with genes in RGD; all entries not matching exactly one active rgd id, are skipped
            Iterator<HOrthologEntry> it = group.getEntries().iterator();
            while( it.hasNext() ) {
                HOrthologEntry entry = it.next();

                List<Gene> genes = dao.getGenesByGeneId(entry.getGeneId());
                if( genes.isEmpty() ) {
                    noRgdId++;
                    it.remove();
                    continue;
                }
                if( genes.size()>1 ) {
                    multiRgdId++;
                    it.remove();
                    continue;
                }
                entry.setGene(genes.get(0));
            }

            group.buildIncomingList("HOMOLOGENE");

            // load data from rgd
            group.inRgdList = dao.getOrthologsForGroupId(group.getHomolegeneId(), "HOMOLOGENE");

            // determine which entries should be inserted, deleted, updated or untouched
            for( Ortholog incoming: group.incomingList ) {

                boolean wasMatch = false;
                Iterator<Ortholog> it2 = group.inRgdList.iterator();
                while( it2.hasNext() ) {
                    Ortholog inRgd = it2.next();
                    if( uniqueKeyMatch(incoming, inRgd) ) {
                        it2.remove();
                        group.matchList.add(inRgd);
                        wasMatch = true;
                        break;
                    }
                }

                if( !wasMatch ) {
                    group.insertList.add(incoming);
                }
            }

            group.deleteList.addAll(group.inRgdList);
        }
    }

    void loadGroups(Collection<HOrthologGroup> groups) throws Exception {

        int insertCount = 0;
        int deleteCount = 0;
        int updateCount = 0;
        int matchCount = 0;

        for( HOrthologGroup group: groups ) {

            if( !group.insertList.isEmpty() ) {
                dao.insertOrthologs(group.insertList);
                insertCount += group.insertList.size();
            }

            if( !group.deleteList.isEmpty() ) {
                dao.deleteOrthologs(group.deleteList);
                deleteCount += group.deleteList.size();
            }

            if( !group.matchList.isEmpty() ) {
                dao.updateLastModified(group.matchList, 70);
                matchCount += group.matchList.size();
            }

            if( !group.updateList.isEmpty() ) {
                throw new Exception("ORTHOLOG updates are not supported");
            }
        }

        System.out.println("INSERTED "+insertCount);
        System.out.println("DELETED "+deleteCount);
        System.out.println("UPDATED "+updateCount);
        System.out.println("MATCHED "+matchCount);
    }


    // must match by SRC_RGD_ID, DEST_RGD_ID
    boolean uniqueKeyMatch(Ortholog o1, Ortholog o2) {
        return o1.getSrcRgdId()==o2.getSrcRgdId() &&
                o1.getDestRgdId()==o2.getDestRgdId();
    }

    class HOrthologEntry {

        private int groupId;
        private int homologeneId;
        private int taxonomicId;
        private String geneId;
        private String dataSet; // source

        private int speciesTypeKey; // derived from taxonomicId
        private Gene gene; // gene matching RGD ID

        public int getHomologeneId() {
            return homologeneId;
        }

        public void setHomologeneId(int homologeneId) {
            this.homologeneId = homologeneId;
        }

        public int getTaxonomicId() {
            return taxonomicId;
        }

        public void setTaxonomicId(int taxonomicId) {
            this.taxonomicId = taxonomicId;

            if( taxonomicId==9606 )
                this.speciesTypeKey = SpeciesType.HUMAN;
            else if( taxonomicId==10090 )
                this.speciesTypeKey = SpeciesType.MOUSE;
            else if( taxonomicId==10116 )
                this.speciesTypeKey = SpeciesType.RAT;
            else
                this.speciesTypeKey = SpeciesType.ALL;
        }

        public String getGeneId() {
            return geneId;
        }

        public void setGeneId(String geneId) {
            this.geneId = geneId;
        }

        public String getDataSet() {
            return dataSet;
        }

        public void setDataSet(String dataSet) {
            this.dataSet = dataSet;
        }

        public Gene getGene() {
            return gene;
        }

        public void setGene(Gene gene) {
            this.gene = gene;
        }

        public int getSpeciesTypeKey() {
            return speciesTypeKey;
        }
    }

    class HOrthologGroup {

        private int homolegeneId;
        private List<HOrthologEntry> entries = new ArrayList<HOrthologEntry>();

        public List<Ortholog> incomingList = new ArrayList<Ortholog>();
        public List<Ortholog> inRgdList;

        public List<Ortholog> insertList = new ArrayList<Ortholog>();
        public List<Ortholog> deleteList = new ArrayList<Ortholog>();
        public List<Ortholog> updateList = new ArrayList<Ortholog>();
        public List<Ortholog> matchList = new ArrayList<Ortholog>();

        // syncer for associations: weak references and paralogs
        private AssociationSyncer assocSyncer = new AssociationSyncer();

        public int getHomolegeneId() {
            return homolegeneId;
        }

        public void setHomolegeneId(int homolegeneId) {
            this.homolegeneId = homolegeneId;
        }

        public List<HOrthologEntry> getEntries() {
            return entries;
        }

        public void setEntries(List<HOrthologEntry> entries) {
            this.entries = entries;
        }

        void buildIncomingList(String xrefDataSource) throws CloneNotSupportedException {

            HOrthologEntry entry1, entry2;
            for( int i=0; i<entries.size()-1; i++ ) {
                entry1 = entries.get(i);

                for( int j=i+1; j<entries.size(); j++ ) {
                    entry2 = entries.get(j);

                    // use longer data set
                    String dataSet = entry1.getDataSet();
                    if( entry2.getDataSet().length() > dataSet.length() )
                        dataSet = entry2.getDataSet();

                    Ortholog o = new Ortholog();
                    o.setXrefDataSet(dataSet);
                    o.setOrthologTypeKey(11);
                    o.setXrefDataSrc(xrefDataSource);
                    o.setCreatedBy(70);
                    o.setLastModifiedBy(70);
                    o.setCreatedDate(new Date());
                    o.setLastModifiedDate(new Date());

                    o.setGroupId(getHomolegeneId());
                    o.setSrcRgdId(entry1.getGene().getRgdId());
                    o.setSrcSpeciesTypeKey(entry1.getSpeciesTypeKey());
                    o.setDestRgdId(entry2.getGene().getRgdId());
                    o.setDestSpeciesTypeKey(entry2.getSpeciesTypeKey());
                    incomingList.add(o);

                    Ortholog o2 = (Ortholog) o.clone();
                    o2.setDestRgdId(entry1.getGene().getRgdId());
                    o2.setDestSpeciesTypeKey(entry1.getSpeciesTypeKey());
                    o2.setSrcRgdId(entry2.getGene().getRgdId());
                    o2.setSrcSpeciesTypeKey(entry2.getSpeciesTypeKey());
                    incomingList.add(o2);
                }
            }
        }

    }
}
