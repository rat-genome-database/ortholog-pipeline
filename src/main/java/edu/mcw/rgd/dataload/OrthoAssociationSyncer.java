package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.sync.AssociationSyncer;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 7/22/13
 * Time: 11:44 AM
 * Synchronizes associations; same as AssociationSyncer from rgdcore, but it does not optimize indels
 */
public class OrthoAssociationSyncer extends AssociationSyncer {

    @Override
    public boolean cloneObjectForIndelOptimization(Object objInsert, Object objDelete) {
        return false;
    }
}
