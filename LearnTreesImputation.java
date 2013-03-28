
import net.maizegenetics.pal.alignment.Alignment;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * Learning algorithms for imputation
 * @author kelly
 */
public class LearnTreesImputation {
    
    public static byte[][] GetAllelesTwoSites(Alignment a, int site1, int site2) { //returns a byte 2-D array of alleles at the input sites
        byte[][] allelesInByte= new byte[2][a.getSequenceCount()];
        for (int i=0; i<a.getSequenceCount();i++) {
            allelesInByte[0][i]= a.getBase(i, site1);
            allelesInByte[1][i]= a.getBase(i, site2);
        }
        
        return allelesInByte;
    }
    
    public static int[] MaskAlleles(byte allele, byte[] taxaAlleles) {//mask a byte array of allele states based on an allele in byte
        int[] mask= new int[taxaAlleles.length];
        for (int i=0; i<taxaAlleles.length; i++){
            if (taxaAlleles[i] == allele) mask[i]= 1;
        }
        
        return mask;
    }
                
    
    public static int FindBestSite(Alignment a, int siteIndex) { //calculates the information gain with respect to the target site for each other site, returns index of site with greatest gain
        int maxInfo;
        int[] infoGain= new int[a.getSiteCount()];
        for (int i= 0;i<a.getSiteCount();i++) {
            byte[][] compareSites= LearnTreesImputation.GetAllelesTwoSites(a,siteIndex,i);
            int[][] alleleFreqAtr= a.getAllelesSortedByFrequency(i, true);
            
            
        }
        return maxInfo;
    }
    
    public static void PruneTree() { //prune based on accuracy in validation set (maybe add complexity penalty?)
        
    }
}
