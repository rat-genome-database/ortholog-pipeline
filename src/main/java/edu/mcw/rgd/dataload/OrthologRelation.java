package edu.mcw.rgd.dataload;

/*
 * @since Mar 22, 2007
 * @author dli
 */
public class OrthologRelation {
    int srcRgdId;
    int destRgdId;
    int srcSpeciesTypeKey;
    int destSpeciesTypeKey;
    String srcOtherId;
    String destOtherId;
    String dataSource;
    String dataSetName;

    public String getDataSetName() {
        return dataSetName;
    }
    
    public void setDataSetName(String dataSetName) {
        this.dataSetName = dataSetName;
    }
    
    public String getDataSource() {
        return dataSource;
    }
    public void setDataSource(String dataSource) {
        this.dataSource = dataSource;
    }
    public int getDestRgdId() {
        return destRgdId;
    }
    public void setDestRgdId(int destRgdId) {
        this.destRgdId = destRgdId;
    }

    public int getSrcRgdId() {
        return srcRgdId;
    }
    public void setSrcRgdId(int srcRgdId) {
        this.srcRgdId = srcRgdId;
    }
    
    public String getDestOtherId() {
        return destOtherId;
    }
    public void setDestOtherId(String destOtherId) {
        this.destOtherId = destOtherId;
    }
    public String getSrcOtherId() {
        return srcOtherId;
    }
    public void setSrcOtherId(String srcOtherId) {
        this.srcOtherId = srcOtherId;
    }
    
    public int getDestSpeciesTypeKey() {
        return destSpeciesTypeKey;
    }
    public void setDestSpeciesTypeKey(int destSpeciesTypeKey) {
        this.destSpeciesTypeKey = destSpeciesTypeKey;
    }
    public int getSrcSpeciesTypeKey() {
        return srcSpeciesTypeKey;
    }
    
    public void setSrcSpeciesTypeKey(int srcSpeciesTypeKey) {
        this.srcSpeciesTypeKey = srcSpeciesTypeKey;
    }  
    
    public String toString() {
        return "srcRgdId="+srcRgdId+" destRgdId="+destRgdId
                +" srcSpecies="+srcSpeciesTypeKey+" destSpecies="+destSpeciesTypeKey
                +" dataSet="+dataSetName;
    }
}
