
import java.util.Arrays;
import java.util.List;
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
    static boolean slidingWindow;
    static int windowSize;
    static int startSite;
    static int endSite;
    static double windowCov; //harmonic mean
    static double windowSegSites;
    static double[] windowSiteFreq;
    static Alignment alignment;
    static Alignment noMissingAlignment;
    static int segSites;
    static double[] noMissingCoverage;
    static double[] unfoldedSiteFreq;
    static int outSequence= -1;
    static int[][] codedAlignment;
    
    public static void SetReference(String file, String refTaxon,String chrNum) {
        alignment= ImportUtils.readFromHapmap(file,chrNum);
        System.out.println("fdsf");
        System.out.println("Number of sites: "+alignment.getSiteCount()+"/n"+"Number of loci: "+alignment.getSiteCount());
        for (int index= 0;index < alignment.getSequenceCount();index++) {
            if (alignment.getFullTaxaName(index).equals(refTaxon)) outSequence= index;
        }
        if (outSequence != -1) {
            System.out.println("Number of sites in alignment: "+alignment.getNumLoci()+"/n"+
            "Number of sequences in alignment: "+alignment.getSequenceCount()+"/n"+
            "Outgroup: "+alignment.getFullTaxaName(outSequence));
        }else{
            System.out.println("reference taxon indicated not found in alignment file");
        }
    }
    
    public static void RemoveSitesMissingInOutgroup() {
        MutableSingleEncodeAlignment mna= new MutableSingleEncodeAlignment(alignment,alignment.getSequenceCount(),alignment.getSiteCount());
        for (int i= 0;i < mna.getNumLoci();i++) {
            if (mna.getBase(outSequence, i) == 'N') mna.removeSite(i);
            }
        mna.clean();
        noMissingAlignment= mna;
        //remove a that are missing in the outgroup
        //generate a new alignment to hold the sequences of interest (minus missing a) but still including the outgroup
        //this is ok because when coded the outgroup will always equal zero
    }
    
    public static void RecodeAlignment() {
        codedAlignment= new int[noMissingAlignment.getSiteCount()][noMissingAlignment.getSequenceCount()];
        for (int i= 0; i < noMissingAlignment.getSiteCount(); i++) {
            double cov= 0;
            if (noMissingAlignment.isPolymorphic(i) == true) segSites++;
            for (int j= 0; j < noMissingAlignment.getSequenceCount(); j++) {
                if (noMissingAlignment.getBase(j, i) == 'N') cov++;
                else if (noMissingAlignment.getBase(outSequence, i) == noMissingAlignment.getBase(j, i)) codedAlignment[i][j]= 0;
                else if (noMissingAlignment.getBase(j, i) == 'A' ||noMissingAlignment.getBase(j, i) == 'T'
                            ||noMissingAlignment.getBase(j, i) == 'C'||noMissingAlignment.getBase(j, i) == 'G') {
                    codedAlignment[i][j]= 2;
                    unfoldedSiteFreq[i]+= 2;
                }
                else {
                    codedAlignment[i][j]= 1;
                    unfoldedSiteFreq[i]++;
                }
            }
            noMissingCoverage[i]= cov/(double) noMissingAlignment.getSequenceCount();
        }
        //recode alignment so that derived variants (different allele than outgroup) equals 2, hets equal 1 and missing and ancestral equals 0
        //also, make a vector that holds the coverage data
        //can use this to calculate all summary statistics
    }
    
    public static void WindowStats() {
        windowCov= 0;
        windowSegSites= 0;
        //find harmonic mean of the coverage for this window
        for (int i= startSite;i < endSite; i++) {
            windowCov+= 1/noMissingCoverage[i];
            if (noMissingAlignment.isPolymorphic(i) == true) windowSegSites++;
        }
        //summarize the unfolded site frequency spectrum for the window
        double[] xi= Arrays.copyOfRange(unfoldedSiteFreq,startSite,endSite);
        Arrays.sort(xi);
        List listXi= Arrays.asList(xi);
        for (int j= 0;j < noMissingAlignment.getSequenceCount();j++) {
            windowSiteFreq[j]= listXi.lastIndexOf(j)-listXi.indexOf(j);
        }
    }
    
    public static double ThetaW() {
        double a= 0;
        //scale 2n (number of chromosomes) by coverage
        for (double n= 1;n < 2*noMissingAlignment.getSequenceCount()*windowCov; n++) {
            a+=1/n;
        }
        return windowSegSites/a;
    }
    
    public void ThetaPi(int[][] codedAlignment) {
        
    }
    
    public void ThetaH(int[][] codedAlignment) {
        
    }
    
    public static void main (String args[]) {
        slidingWindow= true;
        windowSize= 10;
        String fileRoot= "/Users/kelly/Documents/GBS/FinalRev1_BPECFilteredSNPsSubset/";
        SetReference(fileRoot+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr10_NAMTripsacum.hmp.txt.gz","tripsacum:C08L7ACXX:6:250048015","10");
        RemoveSitesMissingInOutgroup();
        RecodeAlignment();
        
    }
    
}
