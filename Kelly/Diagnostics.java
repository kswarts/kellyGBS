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
        long[][] depth= new long[2][200];//first array (index 0) is for masked sites, index 1 for unmasked sites
        int which= 0;
        int currDepth= 0;
        int unmaskedSite;
        for (int site = 0; site < key.getSiteCount(); site++) {
            unmaskedSite= unmasked.getSiteOfPhysicalPosition(key.getPositionInLocus(site), key.getLocus(site));
            which= (key.getMajorAllele(site)==Alignment.UNKNOWN_ALLELE)?1:0;
            for (int taxon = 0; taxon < key.getSequenceCount(); taxon++) {
                currDepth= 0;
                for (byte i:unmasked.getDepthForAlleles(taxon, unmaskedSite)) {currDepth+= i;}
                depth[which][currDepth]++;
            }
        }
        try {
            File outputFile = new File(keyFile.substring(0, keyFile.indexOf("maskKey")) + "SiteDepthDiagnostic.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("MaskedSites\tUnmaskedSites");
            for (int d = 0; d < depth[0].length; d++) {
                outStream.writeBytes("\n");
                outStream.writeLong(depth[0][d]);
                outStream.writeBytes("\t");
                outStream.writeLong(depth[0][d]);
            }
        }
        catch (Exception e) {
        }
    }
    
    public static void main(String[] args) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        
//        String dir= "/Users/kls283/Documents/Imputation/";
        String dir= "/Users/kls283/Desktop/Imputation/";
        String keyFile= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7Small.hmp.h5";
        String unmasked= dir+"AllZeaGBS_v2.7_SeqToGenos_part11.hmp.h5";
        maskedSitesDepth(keyFile, unmasked);
        
        }
}
