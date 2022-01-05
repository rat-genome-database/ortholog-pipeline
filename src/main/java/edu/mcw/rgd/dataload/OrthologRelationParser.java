package edu.mcw.rgd.dataload;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/*
 * Parse the homolog relationship data from file (human-mouse and human-rat relations)
 * <p>
 * format as of sample lines from file 'human_all_hcop_sixteen_column.txt'
 * <pre>
 * 0:speciesid 1:hu-eg_id 2:hu-ens 3:hu hgnc 4:hu-gene-name 5:hu-gene-sym 6:hu-chr 7:hu-assert-ids
10116	2	ENSG00000175899	HGNC:7	alpha-2-macroglobulin	A2M	12p13.31	EOG76T5XW,ENSP00000323929,ENSG00000175899
 * 8:ortho-eg-id 9:o ens 10:ext id 11:gene name 12:gene sym 13:o-chr 14:o-assert-ids 15:o-assert-sources
297568	ENSRNOG00000042228	RGD:1584999	alpha-1-inhibitor III	A1i3	4	EOG76T5XW,ENSRNOP00000008574,ENSRNOG00000042228	OrthoDB,OrthoMCL,Ensembl
 * </pre>
 */
public class OrthologRelationParser {

    private List<OrthologRelation> orthologRelations;
    int totalRelations=0;
    int skippedRelations=0;
    
    protected final Logger processLog = LogManager.getLogger("process");
    private OrthologRelationFile file;

    public List<OrthologRelation> parse(int speciesTypeKey) throws Exception {
        file.downloadFile(speciesTypeKey);
        readFile(speciesTypeKey);
        List<OrthologRelation> orthoRelations = new ArrayList<>(orthologRelations);

        if( file.isFileFromHCOP() ) {
            file.downloadFileFromNcbi();
            readFile(speciesTypeKey);
            orthoRelations.addAll(orthologRelations);
        }

        processLog.info("  total relations parsed:"+(totalRelations+skippedRelations));
        processLog.info("  human/"+SpeciesType.getCommonName(speciesTypeKey)+" relations: "+totalRelations);
        //processLog.info("  other species relations:         "+skippedRelations);

        if( totalRelations<5000 ) {
            throw new Exception("POSSIBLE PROBLEM WITH SOURCE FILE: only "+totalRelations+" found in the file");
        }
        return orthoRelations;
    }

    /**
     * read file and load orthologs between human and given species
     * @throws Exception
     */
    public void readFile(int speciesTypeKey2) throws Exception {

        BufferedReader br = openInputFile();
        if( br==null ) {
            return;
        }

        String taxId = Integer.toString(SpeciesType.getTaxonomicId(speciesTypeKey2));

        int speciesTypeKey1 = SpeciesType.HUMAN;
        String humanTaxId = Integer.toString(SpeciesType.getTaxonomicId(speciesTypeKey1));

        String line;
        while ((line=br.readLine())!=null) {

            String[] cols = line.split("[\t]", -1);

            if( file.isFileFromHCOP() ) {
                // first line contains species of interest
                String speciesTaxId = cols[0];
                if (!speciesTaxId.equals(taxId)) {
                    skippedRelations++;
                    continue;
                }
                totalRelations++;

                addOrthologFromHcop(cols, speciesTypeKey1, speciesTypeKey2);
            } else {
                // NCBI file line format:
                // #tax_id	GeneID	relationship	Other_tax_id	Other_GeneID
                //
                // first line contains human
                if (!humanTaxId.equals(cols[0])) {
                    skippedRelations++;
                    continue;
                }
                String speciesTaxId = cols[3];
                if (!speciesTaxId.equals(taxId)) {
                    skippedRelations++;
                    continue;
                }
                totalRelations++;

                addOrthologFromNcbi(cols, speciesTypeKey1, speciesTypeKey2);
            }
        }
        br.close();
    }

    BufferedReader openInputFile() throws IOException {
        //System.out.println(getClass());
        orthologRelations=new ArrayList<> ();

        //check if the file exist
        String localFileName = getFile().getLocalFileName();
        File file=new File(localFileName);
        if (!file.exists()) {
            processLog.error("The file "+ localFileName + " doesn't exist");
            return null;
        }

        processLog.info("Parsing file: " + localFileName);

        return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(localFileName))));
    }

    // column with sources for ortholog:
    // OrthoDB,OrthoMCL,Ensembl,Ensembl
    // remove duplicates and sort alphabetically
    String sanitizeXRefDataSet(String xrefDataSet) {

        // compute non-duplicate list of sources from the current XREF_DATA_SET
        String[] sourcesOrig = xrefDataSet.split(",");
        if( sourcesOrig.length>1 ) {
            Set<String> sourcesUnique = new TreeSet<>(Arrays.asList(sourcesOrig));
            return Utils.concatenate(sourcesUnique, ", ");
        } else {
            return sourcesOrig[0];
        }
    }

    void addOrthologFromHcop(String[] cols, int speciesTypeKey1, int speciesTypeKey2) {
        String geneIdHuman = cols[1], geneIdOther = cols[8];

        // column with sources for ortholog:
        // OrthoDB,OrthoMCL,Ensembl,Ensembl
        String xrefDataSet = sanitizeXRefDataSet(cols[15]);

        OrthologRelation entry = new OrthologRelation();
        entry.setDataSource(OrthologRelationLoadingManager.getInstance().getXrefDataSrc());
        entry.setDataSetName(xrefDataSet);
        entry.setSrcSpeciesTypeKey(speciesTypeKey1);
        entry.setDestSpeciesTypeKey(speciesTypeKey2);
        entry.setSrcOtherId(geneIdHuman);
        entry.setDestOtherId(geneIdOther);
        orthologRelations.add(entry);
    }

    void addOrthologFromNcbi(String[] cols, int speciesTypeKey1, int speciesTypeKey2) {
        String geneIdHuman = cols[1], geneIdOther = cols[4];

        OrthologRelation entry = new OrthologRelation();
        entry.setDataSource("NCBI");
        entry.setDataSetName(cols[2]);
        entry.setSrcSpeciesTypeKey(speciesTypeKey1);
        entry.setDestSpeciesTypeKey(speciesTypeKey2);
        entry.setSrcOtherId(geneIdHuman);
        entry.setDestOtherId(geneIdOther);
        orthologRelations.add(entry);
    }

    public int getSkippedRelations() {
        return skippedRelations;
    }
    public int getTotalRelations() {
        return totalRelations;
    }

    public void setFile(OrthologRelationFile file) {
        this.file = file;
    }

    public OrthologRelationFile getFile() {
        return file;
    }
}
