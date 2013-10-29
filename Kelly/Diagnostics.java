/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.prefs.TasselPrefs;

/**
 *
 * @author kls283
 */
//gets site depth statistics from the unmasked file. can restrict sites calculated for by the key file, but takes taxon from the unmasked file
public class Diagnostics {
    public static void maskedSitesDepth(String keyFile, String unmaskedFile) {
        MutableNucleotideAlignmentHDF5 key= MutableNucleotideAlignmentHDF5.getInstance(keyFile);
        MutableNucleotideAlignmentHDF5 unmasked= MutableNucleotideAlignmentHDF5.getInstance(unmaskedFile);
        int[][] depth= new int[2][200];//first array (index 0) is for masked sites, index 1 for unmasked sites
        int which= 0;
        int currDepth= 0;
        int unmaskedSite;
        for (int site = 0; site < key.getSiteCount(); site++) {
            if (depth[0][0]==Long.MAX_VALUE||depth[1][0]==Long.MAX_VALUE) {System.out.println("Reached long max at site "+(site-1)); break;}
            unmaskedSite= unmasked.getSiteOfPhysicalPosition(key.getPositionInLocus(site), key.getLocus(site));
            which= (key.getMajorAllele(site)==Alignment.UNKNOWN_ALLELE)?1:0;
            for (int taxon = 0; taxon < unmasked.getSequenceCount(); taxon++) {
                currDepth= 0;
                for (byte i:unmasked.getDepthForAlleles(taxon, unmaskedSite)) {currDepth+= i;}
                depth[which][currDepth]++;
            }
        }
        System.out.println("Depth\tMaskedSites\tUnmaskedSites");
        for (int i = 0; i < depth[0].length; i++) {
            System.out.println(i+"\t"+depth[0]+"\t"+depth[1]);
        }
        try {
            File outputFile = new File(keyFile.substring(0, keyFile.indexOf(".hmp")) + "SiteDepthDiagnostic.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("Depth\tMaskedSites\tUnmaskedSites");
            for (int d = 0; d < depth[0].length; d++) {
                outStream.writeBytes("\n"+d+"\t");
                outStream.writeLong(depth[0][d]);
                outStream.writeBytes("\t");
                outStream.writeLong(depth[1][d]);
            }
        }
        catch (Exception e) {
        }
    }
    
    public static void main(String[] args) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        
        String dir= "/home/kls283/Documents/Imputation/";
//        String dir= "/Users/kls283/Desktop/Imputation/";
        String keyFile= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7.hmp.h5";
        String unmasked= dir+"AllZeaGBS_v2.7wDepth.hmp.h5";
        maskedSitesDepth(keyFile, unmasked);
        
        }
}