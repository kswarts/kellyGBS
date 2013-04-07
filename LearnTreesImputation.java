
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import net.maizegenetics.pal.alignment.*;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * Learning algorithms for imputation
 * @author kelly
 */
public class LearnTreesImputation {
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte('N');
    public static byte[] genotypicN= new byte[]{diploidN,diploidN};
    
    public static String[] GetDiploidAllelesAtSite(Alignment a, int site) {
        String[] alleles= new String[a.getSequenceCount()];
        for (int i=0; i<a.getSequenceCount();i++) {
            alleles[i]= a.getBaseAsString(i, site);
        }
        return alleles;
    }
    
        public static String[] GetHaploidAllelesAtSite(Alignment a, int site) {
        String[] alleles= new String[a.getSequenceCount()];
        byte[]hap;
        for (int i=0; i<a.getSequenceCount();i++) {
            hap= a.getDepthForAllele(i, site);
            alleles[i]= NucleotideAlignmentConstants.getNucleotideIUPAC(hap[0])+NucleotideAlignmentConstants.getNucleotideIUPAC(hap[1]);
        }
        return alleles;
    }
    
    public static byte[][] GetAllelesTwoSites(Alignment a, int site1, int site2) { //returns a byte 2-D array of alleles at the input sites
        byte[][] allelesInByte= new byte[2][a.getSequenceCount()];
        for (int i=0; i<a.getSequenceCount();i++) {
            allelesInByte[0][i]= a.getBase(i, site1);
            allelesInByte[1][i]= a.getBase(i, site2);
        }
        
        return allelesInByte;
    }
    
    public static int[] MaskAlleles(Alignment a, int site, byte allele) {//mask a site based on an allele in byte (diploid value)
        int[] mask= new int[a.getSequenceCount()];
        for (int i=0; i<mask.length; i++){
            if (a.getBase(i, site) == allele) mask[i]= 1;
        }
        
        return mask;
    }
    
    public static Alignment SubsetAlignmentByAllele(Alignment a, int site, byte allele) { //make a subset alignment that only contains sequences where a site specified has a specific diploid allele value
        MutableNucleotideAlignment align= MutableNucleotideAlignment.getInstance(a);
        int[] mask= LearnTreesImputation.MaskAlleles(a,site,allele);
        for (int taxa= 0; taxa<a.getSequenceCount(); taxa++) {
            if (mask[taxa]==0) align.clearSiteForRemoval(site);
        }
        align.clean();
        
        return align;
    }
    
    public static Alignment SelectPresent(Alignment a, int site) {
        MutableNucleotideAlignment align= MutableNucleotideAlignment.getInstance(a);
        for (int taxa= 0; taxa<a.getSequenceCount(); taxa++) {
            if (a.getBase(taxa, site) == diploidN) align.clearSiteForRemoval(site);
        }
        align.clean();
        
        return align;
    }
    
    public static double GetEntropy(Alignment a, int siteIndex) {
        double Pmaj= a.getMajorAlleleFrequency(siteIndex);
        double Pmin= a.getMinorAlleleFrequency(siteIndex);
        double Pmiss= 1-((double) a.getTotalGametesNotMissing(siteIndex)/((double) a.getSequenceCount()*2));
        double H= (a.getMinorAllele(siteIndex)==diploidN?0:Pmaj*-Math.log10(Pmaj))+(Pmin*-Math.log10(Pmin))+(Pmiss*-Math.log10(Pmiss));
        
        return H;
    }
    
    public static double GetConditionalEntropy(Alignment a, int site1, int site2, byte allele) { //site one is the focal site which allele states are conditioned on
        Alignment align= LearnTreesImputation.SubsetAlignmentByAllele(a, site1, allele);
        double condH= LearnTreesImputation.GetEntropy(align,site2);
        
        return condH;
    }
                
    
    public static int FindBestSite(Alignment a, int siteIndex) { //calculates the information gain with respect to the target site for each other site, returns index of site with greatest gain        
        double maxInfo= 0;
        int bestSiteIndex= siteIndex;
        double Pmaj= a.getMajorAlleleFrequency(siteIndex);
        byte majAllele= a.getMajorAllele(siteIndex);
        double Pmin= a.getMinorAlleleFrequency(siteIndex);
        byte minAllele= a.getMinorAllele(siteIndex);
        double[] infoGain= new double[a.getSiteCount()];
//        System.out.println("index site: "+a.getSNPID(siteIndex)+"\n");
        
        for (int site2= 0;site2<a.getSiteCount();site2++) {
            infoGain[site2]= Pmaj*LearnTreesImputation.GetConditionalEntropy(a,siteIndex,site2,majAllele)+Pmin*LearnTreesImputation.GetConditionalEntropy(a,siteIndex,site2,minAllele);
            if (infoGain[site2]>maxInfo) bestSiteIndex= site2;
//            System.out.println(infoGain[site2]);
        }
        try{
            String outFile= "/Users/kelly/Documents/GBS/FinalRev1_BPECFilteredSNPsSubset/FindBestSiteTest.txt";
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 655360));
          outStream.writeBytes("SNPID\tInfoGain\tSiteCov\tMAF");
          
          for (int site= 0; site < a.getSiteCount(); site++) {              
              outStream.writeBytes("\n"+a.getSNPID(site)+"\t"+infoGain[site]+"\t"+(double)a.getTotalGametesNotMissing(site)/((double)a.getSequenceCount()*2)+"\t"+a.getMajorAlleleFrequency(site));              
           }
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
        
        return bestSiteIndex;
    }
    
    public static void PruneTree() { //prune based on accuracy in validation set (maybe add complexity penalty?)
        
    }
    
    public static void main (String args[]) {
       
       String dir= "/Users/kelly/Documents/GBS/FinalRev1_BPECFilteredSNPsSubset/";
       String inFile= "AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr10_NAM3Polymorphic";
       Alignment a= ImportUtils.readFromHapmap(dir+inFile+".hmp.txt", null);
       int refSite= 1;
       String SNP= a.getSNPID(refSite);
       String[] alleles= LearnTreesImputation.GetHaploidAllelesAtSite(a, refSite);
       Alignment noMissing= LearnTreesImputation.SelectPresent(a, refSite);
       double H= LearnTreesImputation.GetEntropy(noMissing, refSite);
       System.out.println("H for "+SNP+": "+H);
       for (int i=0;i<alleles.length;i++){System.out.println(alleles[i]);}
//       int root= LearnTreesImputation.FindBestSite(noMissing, refSite);
//       System.out.println("Best Site: "+a.getSNPID(root)+" ("+root+")");
       


   }
}
