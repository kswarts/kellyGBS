/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.prefs.TasselPrefs;
import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author kls283
 */
//gets site depth statistics from the unmasked file. can restrict sites calculated for by the key file, but takes taxon from the unmasked file
public class Diagnostics {
    public static void maskedSitesDepth(String keyFile, String unmaskedFile) {
        MutableNucleotideAlignmentHDF5 key= MutableNucleotideAlignmentHDF5.getInstance(keyFile);
        MutableNucleotideAlignmentHDF5 unmasked= MutableNucleotideAlignmentHDF5.getInstance(unmaskedFile);
        long[][] depth= new long[2][1000];//first array (index 0) is for masked sites, index 1 for unmasked sites
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
    
    public static void removeIndelsForBeagle(String[] inFiles) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        for (String file:inFiles) {
            Alignment a= ImportUtils.readGuessFormat(file, true);
            ArrayList<Integer> keepSites= new ArrayList<>();
            for (int site = 0; site < a.getSiteCount(); site++) {
                if (a.getMajorAllele(site)!= NucleotideAlignmentConstants.GAP_ALLELE&&
                        a.getMajorAllele(site)!= NucleotideAlignmentConstants.INSERT_ALLELE&&
                        a.getMinorAllele(site)!= NucleotideAlignmentConstants.GAP_ALLELE&&
                        a.getMinorAllele(site)!= NucleotideAlignmentConstants.INSERT_ALLELE)
                    keepSites.add(site);
            }
            FilterAlignment fa= FilterAlignment.getInstance(a, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));
            ExportUtils.writeToVCF(fa, file.substring(0, file.indexOf(".vcf"))+"NoIndels.vcf.gz", '\t');
        }
    }
    
    public static void removeIndelsForBeagle(String dir) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        File[] inFiles= new File(dir).listFiles();
        
        for (File file:inFiles) {
            String readFile= null;
            if (file.isFile() && file.getName().endsWith(".vcf.gz")) readFile= file.getName();
            if (readFile==null) continue;
            Alignment a= ImportUtils.readGuessFormat(readFile, true);
            ArrayList<Integer> keepSites= new ArrayList<>();
            for (int site = 0; site < a.getSiteCount(); site++) {
                if (a.getMajorAllele(site)!= NucleotideAlignmentConstants.GAP_ALLELE&&
                        a.getMajorAllele(site)!= NucleotideAlignmentConstants.INSERT_ALLELE&&
                        a.getMinorAllele(site)!= NucleotideAlignmentConstants.GAP_ALLELE&&
                        a.getMinorAllele(site)!= NucleotideAlignmentConstants.INSERT_ALLELE)
                    keepSites.add(site);
            }
            FilterAlignment fa= FilterAlignment.getInstance(a, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));
            ExportUtils.writeToVCF(fa, readFile.substring(0, readFile.indexOf(".vcf"))+"NoIndels.vcf.gz", '\t');
        }
    }
    
    public static void main(String[] args) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        String dir;
        
        dir= "/home/kls283/Documents/Imputation/";
//        String dir= "/Users/kls283/Desktop/Imputation/";
        String keyFile= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7.hmp.h5";
        String unmasked= dir+"AllZeaGBS_v2.7wDepth.hmp.h5";
        maskedSitesDepth(keyFile, unmasked);
        
        dir= "/home/kls283/Documents/Imputation/";
        String[] files= new String[] {dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy282Allchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy282Allchr8.vcf.gz"};
        removeIndelsForBeagle(dir);
        removeIndelsForBeagle(files);
        
        }
}