
import net.maizegenetics.pal.alignment.*;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author kelly
 */
public class TestsOfSelection {
    Alignment alignment;
    Alignment noMissingAlignment;
    int outSequence;
    int[] coverage;
    int[][] codedAlignment;
    
    public void setReference(String file, String refTaxon) {
        try {
            this.alignment= ReadSequenceAlignmentUtils.readBasicAlignments(file, 30);
        } catch (Exception e) {
            //do nothing
        }
        int index= 0;
        while (alignment.getFullTaxaName(index) != refTaxon && index < alignment.getSequenceCount()) {
            index++;
        }
        if (alignment.getFullTaxaName(index) == refTaxon) {
            outSequence= index;
            System.out.println("Number of sites in alignment: "+alignment.getNumLoci()+"/n"+
                    "Number of sequences in alignment: "+alignment.getSequenceCount()+"/n"+
                    "Outgroup: "+alignment.getFullTaxaName(outSequence));
        }else{
            System.out.println("reference taxon indicated not found in alignment file");
    } 
    }
    
    public void removeSitesMissingInOutgroup(Alignment a, int outSequence) {
        this.alignment= a;
        this.outSequence= outSequence;
        //remove sites that are missing in the outgroup
        //generate a new alignment to hold the sequences of interest (minus missing sites) but still including the outgroup
        //this is ok because when coded the outgroup will always equal zero
    }
    
    public void recodeAlignment() {
        //recode alignment so that derived variants (different allele than outgroup) equals 2, hets equal 1 and missing and ancestral equals 0
        //also, make a vector that holds the coverage data
        //can use this to calculate all summary statistics
    }
    
}
