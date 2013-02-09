
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
    static Alignment alignment;
    static Alignment noMissingAlignment;
    static int outSequence= -1;
    static int[] coverage;
    static int[][] codedAlignment;
    
    public static void SetReference(String file, String refTaxon) {
        try {
            alignment= ReadSequenceAlignmentUtils.readBasicAlignments(file, 30);
        } catch (Exception e) {
            //do nothing
        }
        for (int index= 0;index <= alignment.getSequenceCount();index++) {
            if (alignment.getFullTaxaName(index) == refTaxon) outSequence= index;
        }
        if (outSequence != -1) {
            System.out.println("Number of sites in alignment: "+alignment.getNumLoci()+"/n"+
            "Number of sequences in alignment: "+alignment.getSequenceCount()+"/n"+
            "Outgroup: "+alignment.getFullTaxaName(outSequence));
        }else{
            System.out.println("reference taxon indicated not found in alignment file");
    }
    
    public static void RemoveSitesMissingInOutgroup(Alignment a, int outgroup) {
        alignment= a;
        outSequence= outgroup;
        //remove sites that are missing in the outgroup
        //generate a new alignment to hold the sequences of interest (minus missing sites) but still including the outgroup
        //this is ok because when coded the outgroup will always equal zero
    }
    
    public void RecodeAlignment() {
        //recode alignment so that derived variants (different allele than outgroup) equals 2, hets equal 1 and missing and ancestral equals 0
        //also, make a vector that holds the coverage data
        //can use this to calculate all summary statistics
    }
    
    public void ThetaW(int[][] codedAlignment) {
        
    }
    
    public void ThetaPi(int[][] codedAlignment) {
        
    }
    
    public void ThetaH(int[][] codedAlignment) {
        
    }
    
    public static void main (String args[]) {
        String fileRoot= "/Users/kelly/Documents/GBS/FinalRev1_BPECFilteredSNPsSubset/";
        SetReference(fileRoot+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr10_NAMTripsacum.hmp.txt.gz","tripsacum:C08L7ACXX:6:250048015");
        
    }
    
}
