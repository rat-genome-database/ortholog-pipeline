package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.OrthologDAO;
import edu.mcw.rgd.dao.spring.GeneQuery;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.Ortholog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

public class OrthoTool {

    public static void main( String[] args) throws Exception {

        OrthologDAO odao = new OrthologDAO();
        GeneDAO gdao = new GeneDAO();

        String inputFileName = "h:/NCBI_rat_genes.txt";
        String outputFileName = "h:/NCBI_rat_genes3.txt";
        BufferedReader in = Utils.openReader(inputFileName);
        BufferedWriter out = Utils.openWriter(outputFileName);

        String line;

        out.write("Rat RGDID\tRat GENEID\tCurrent symbol\tCurrent name\t\tMouse RGDID\tMouse symbol\tMouse name\t\tHuman RGDID\tHuman symbol\tHuman name\n");

        int mouseOrthologs = 0;
        int humanOrthologs = 0;
        int i=0;
        while( (line=in.readLine())!=null ) {

            String[] cols = line.split("[\\t]", -1);

            int ratRgdId = 0;
            try {
                ratRgdId = Integer.parseInt(getText(cols[0]));
            } catch( NumberFormatException e ) {
                continue;
            }
            int ratGeneId = Integer.parseInt(getText(cols[1]));
            String ratGeneSymbol = getText(cols[2]);
            String ratGeneName = getText(cols[3]);

            String ratData = ratRgdId+"\t"+ratGeneId+"\t"+ratGeneSymbol+"\t"+ratGeneName+"\t\t";

            // mouse ortholog
            List<Gene> mouseGenes;
            List<Gene> humanGenes;
            if( false ) {
                List<Ortholog> orthos = odao.getOrthologsForSourceRgdId(ratGeneId);
                mouseGenes = new ArrayList<>();
                for (Ortholog o : orthos) {
                    if (o.getDestSpeciesTypeKey() == SpeciesType.MOUSE) {
                        mouseGenes.add(gdao.getGene(o.getDestRgdId()));
                    }
                }
                humanGenes = new ArrayList<>();
                for (Ortholog o : orthos) {
                    if (o.getDestSpeciesTypeKey() == SpeciesType.HUMAN) {
                        humanGenes.add(gdao.getGene(o.getDestRgdId()));
                    }
                }
            } else {
                mouseGenes = getAgrMouseGene(ratRgdId, gdao);
                humanGenes = getAgrHumanGene(ratRgdId, gdao);
            }

            do {

                Gene mouseGene = null;
                if( !mouseGenes.isEmpty() ) {
                    mouseGene = mouseGenes.remove(0);
                }
                out.write(ratData);
                if( mouseGene==null ) {
                    out.write("\t\t\t\t");
                } else {
                    mouseOrthologs++;
                    out.write(Integer.toString(mouseGene.getRgdId()));
                    out.write("\t");
                    out.write(mouseGene.getSymbol());
                    out.write("\t");
                    out.write(mouseGene.getName());
                    out.write("\t\t");
                }


                // human ortholog
                Gene humanGene = null;
                if( !humanGenes.isEmpty() ) {
                    humanGene = humanGenes.remove(0);
                }
                if( humanGene==null ) {
                    out.write("\t\t\n");
                } else {
                    humanOrthologs++;
                    out.write(Integer.toString(humanGene.getRgdId()));
                    out.write("\t");
                    out.write(humanGene.getSymbol());
                    out.write("\t");
                    out.write(humanGene.getName());
                    out.write("\n");
                }
            } while( mouseGenes.size() + humanGenes.size() > 0 );

            System.out.println(++i + ".");
        }

        in.close();
        out.close();
        System.out.println("mouse orthos: "+mouseOrthologs+", human orthos: "+humanOrthologs);
    }

    static String getText(String col) {
        if( col.startsWith("\"") ) {
            return col.substring(1, col.length()-1).trim();
        }
        return col.trim();
    }

    static List<Gene> getAgrMouseGene(int ratGeneRgdId, GeneDAO gdao) throws Exception {
        return getAgrGenes(ratGeneRgdId, gdao, 2);
    }

    static List<Gene> getAgrHumanGene(int ratGeneRgdId, GeneDAO gdao) throws Exception {
        return getAgrGenes(ratGeneRgdId, gdao, 1);
    }

    static List<Gene> getAgrGenes(int ratGeneRgdId, GeneDAO gdao, int speciesTypeKey) throws Exception {
        String sql = "SELECT g.*,i.species_type_key FROM agr_orthologs o,genes g,rgd_ids i WHERE gene_rgd_id_1=? "+
                "AND gene_rgd_id_2=g.rgd_id AND g.rgd_id=i.rgd_id AND i.species_type_key=?";
        List<Gene> genes = GeneQuery.execute(gdao, sql, ratGeneRgdId, speciesTypeKey);
        if( genes.isEmpty() ) {
            return genes;
        }

        List<String> isBestScore = new ArrayList<>();
        List<String> isBestRevScore = new ArrayList<>();
        for( int i=0; i<genes.size(); i++ ) {
            Gene g = genes.get(i);
            isBestScore.add(gdao.getStringResult("SELECT is_best_score FROM agr_orthologs WHERE gene_rgd_id_1=? AND gene_rgd_id_2=?", ratGeneRgdId, g.getRgdId()));
            isBestRevScore.add(gdao.getStringResult("SELECT is_best_rev_score FROM agr_orthologs WHERE gene_rgd_id_1=? AND gene_rgd_id_2=?", ratGeneRgdId, g.getRgdId()));
        }

        // count how many of the results have is_best_score='Y' and is_best_rev_score='Y'
        int count = 0;
        for( int i=0; i<genes.size(); i++ ) {
            if( isBestScore.get(i).equals("Y") && isBestRevScore.get(i).equals("Y") ) {
                count++;
            }
        }
        if( count>0 ) {
            // remove non-best scores from the results
            for( int i=genes.size()-1; i>=0; i-- ) {
                if( !(isBestScore.get(i).equals("Y") && isBestRevScore.get(i).equals("Y")) ) {
                    genes.remove(i);
                }
            }
        }
        return genes;
    }
}
