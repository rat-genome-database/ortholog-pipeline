package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Ortholog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;

import java.util.*;

/**
 * @author mtutaj
 * @since 5/14/13
 */
public class OrthologGroup {

    // all related relationships
    public List<OrthologRelation> relations = new ArrayList<>();

    public List<Ortholog> incomingList = new ArrayList<>();

    public boolean add(OrthologRelation rel) {
        if( relations.isEmpty() ) {
            return relations.add(rel);
        }

        // merge duplicates
        for( OrthologRelation r: relations ) {
            if( r.getSrcRgdId()==rel.getSrcRgdId() && r.getDestRgdId()==rel.getDestRgdId() ) {
                // if dataSource is the same, merge dataSetNames
                if( r.getDataSource().equals(rel.getDataSource()) ) {
                    String newDataSetName = mergeDataSetNames(r.getDataSetName(), rel.getDataSetName());
                    if( !newDataSetName.equals(r.getDataSetName()) ) {
                        r.setDataSetName(newDataSetName);
                    }
                    return false;
                } else {
                    // different data sources: current code designed only to merge with NCBI
                    if( !rel.getDataSource().equals("NCBI") ) {
                        throw new RuntimeException("OrthologGroup: cannot merge relations");
                    } else {
                        // merge NCBI with HGNC
                        String newDataSetName = mergeDataSetNames(r.getDataSetName(), rel.getDataSource());
                        if( !newDataSetName.equals(r.getDataSetName()) ) {
                            r.setDataSetName(newDataSetName);
                        }
                        return false;
                    }
                }
            }
        }
        return relations.add(rel);
    }

    String mergeDataSetNames(String dsn1, String dsn2 ) {
        Set<String> set = new TreeSet<>();
        String[] ds = dsn1.split("\\,\\s");
        set.addAll(Arrays.asList(ds));
        ds = dsn2.split("\\,\\s");
        set.addAll(Arrays.asList(ds));
        return Utils.concatenate(set, ", ");
    }

    /**
     * all relations are supposed to be human to other species
     * this method 'complements' the relations by creating all cross relations between species
     */
    public void buildComplementaryRelations() throws Exception {

        Set<OrthologRelation> assocSet = new HashSet<>();

        // create complementary relations
        for( OrthologRelation relHuman: this.relations ) {

            if( relHuman.getSrcSpeciesTypeKey()!=SpeciesType.HUMAN ) {
                throw new Exception("unexpected species type key");
            }

            // create reverse relation
            OrthologRelation rel = new OrthologRelation();
            rel.setSrcRgdId(relHuman.getDestRgdId());
            rel.setDestRgdId(relHuman.getSrcRgdId());
            rel.setSrcSpeciesTypeKey(relHuman.getDestSpeciesTypeKey());
            rel.setDestSpeciesTypeKey(relHuman.getSrcSpeciesTypeKey());
            rel.setSrcOtherId(relHuman.getDestOtherId());
            rel.setDestOtherId(relHuman.getSrcOtherId());
            rel.setDataSource(relHuman.getDataSource());
            rel.setDataSetName(relHuman.getDataSetName());
            assocSet.add(rel);
        }

        relations.addAll(assocSet);
    }
}
