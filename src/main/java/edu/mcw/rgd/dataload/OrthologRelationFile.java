package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.FileDownloader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Set;

/*
 * download a file from external ftp site and store a local copy of that file
 * <p> depending on the species, the source file is either downloaded from HGNC HCOP, or from NCBI
 */
public class OrthologRelationFile extends FileDownloader {

    protected final Logger processLog = LogManager.getLogger("process");

    private boolean isFileFromHCOP;
    private String localFileName;
    private String externalHcopFile;
    private String externalNcbiFile;
    private String localHcopFile;
    private String localNcbiFile;
    private Set<String> hcopSpecies;

    // download a file
    public void downloadFile(int speciesTypeKey) throws Exception {

        // set external and local files depending on species
        String species = SpeciesType.getCommonName(speciesTypeKey);
        if( getHcopSpecies().contains(species) ) {
            setExternalFile(getExternalHcopFile());
            setLocalFile(getLocalHcopFile());
            setFileFromHCOP(true);
            processLog.info("Downloading file " + getExternalFile());
            setLocalFileName(downloadNew());
            processLog.debug("File has been downloaded to "+getLocalFileName());

        } else {
            downloadFileFromNcbi();
        }
    }

    public void downloadFileFromNcbi() throws Exception {

        setExternalFile(getExternalNcbiFile());
        setLocalFile(getLocalNcbiFile());
        setFileFromHCOP(false);

        processLog.info("Downloading file " + getExternalFile());
        setLocalFileName(downloadNew());
        processLog.debug("File has been downloaded to "+getLocalFileName());
    }

    public String getLocalFileName() {
        return localFileName;
    }

    public void setLocalFileName(String localFileName) {
        this.localFileName = localFileName;
    }

    public void setExternalHcopFile(String externalHcopFile) {
        this.externalHcopFile = externalHcopFile;
    }

    public String getExternalHcopFile() {
        return externalHcopFile;
    }

    public void setExternalNcbiFile(String externalNcbiFile) {
        this.externalNcbiFile = externalNcbiFile;
    }

    public String getExternalNcbiFile() {
        return externalNcbiFile;
    }

    public void setLocalHcopFile(String localHcopFile) {
        this.localHcopFile = localHcopFile;
    }

    public String getLocalHcopFile() {
        return localHcopFile;
    }

    public void setLocalNcbiFile(String localNcbiFile) {
        this.localNcbiFile = localNcbiFile;
    }

    public String getLocalNcbiFile() {
        return localNcbiFile;
    }

    public boolean isFileFromHCOP() {
        return isFileFromHCOP;
    }

    public void setFileFromHCOP(boolean fileFromHCOP) {
        isFileFromHCOP = fileFromHCOP;
    }

    public void setHcopSpecies(Set<String> hcopSpecies) {
        this.hcopSpecies = hcopSpecies;
    }

    public Set<String> getHcopSpecies() {
        return hcopSpecies;
    }
}
