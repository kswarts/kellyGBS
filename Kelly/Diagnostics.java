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
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
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
        boolean reachEnd= false;
        MutableNucleotideAlignmentHDF5 key= MutableNucleotideAlignmentHDF5.getInstance(keyFile);
        MutableNucleotideAlignmentHDF5 unmasked= MutableNucleotideAlignmentHDF5.getInstance(unmaskedFile);
        int[][] depth= new int[2][1000];//first array (index 0) is for masked sites, index 1 for unmasked sites
        int which= 0;
        int currDepth= 0;
        int unmaskedSite;
        int inc= Character.getNumericValue(keyFile.charAt(keyFile.indexOf(".hmp")-1));
        for (int site = 0; site < key.getSiteCount(); site+= 10) {
            unmaskedSite= unmasked.getSiteOfPhysicalPosition(key.getPositionInLocus(site), key.getLocus(site));
            which= (key.getMajorAllele(site)==Alignment.UNKNOWN_ALLELE)?1:0;
            for (int taxon = 0; taxon < unmasked.getSequenceCount(); taxon++) {
                currDepth= 0;
                for (byte i:unmasked.getDepthForAlleles(taxon, unmaskedSite)) {currDepth+= i;}
                if (depth[which][currDepth]==Integer.MAX_VALUE) {System.out.println("Reached long max at site "+(site-1)); reachEnd= true; break;}
                else depth[which][currDepth]++;
            }
            if (reachEnd) break;
            System.out.println("Complete site "+site+" of "+key.getSiteCount());
        }
        System.out.println("Depth\tMaskedSites\tUnmaskedSites");
        for (int i = 0; i < depth[0].length; i++) {
            System.out.println(i+"\t"+depth[0][i]+"\t"+depth[1][i]);
        }
        try {
            File outputFile = new File(keyFile.substring(0, keyFile.indexOf(".hmp")) + "SiteDepthDiagnostic.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("Depth\tMaskedSites\tUnmaskedSites");
            for (int d = 0; d < depth[0].length; d++) {outStream.writeBytes("\n"+d+"\t"+depth[0][d]+"\t"+depth[1][d]);}
        }
        catch (Exception e) {
        }
    }
    
    public static void removeIndelsForBeagle(String[] inFiles) {
        for (String file:inFiles) {
            Alignment a= ImportUtils.readGuessFormat(file, true);
            MutableNucleotideAlignment mna= (a.retainsRareAlleles()==true)?MutableNucleotideAlignment.getInstance(a):null;
            ArrayList<Integer> keepSites= new ArrayList<>();
            for (int site = 0; site < a.getSiteCount(); site++) {
                if (a.getMajorAllele(site)== NucleotideAlignmentConstants.GAP_ALLELE||
                        a.getMajorAllele(site)== NucleotideAlignmentConstants.INSERT_ALLELE||
                        a.getMinorAllele(site)== NucleotideAlignmentConstants.GAP_ALLELE||
                        a.getMinorAllele(site)== NucleotideAlignmentConstants.INSERT_ALLELE)
                    continue;
                keepSites.add(site);
                if (a.retainsRareAlleles()&&a.getAlleles(site).length>2) {
                    byte badGeno= AlignmentUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, NucleotideAlignmentConstants.INSERT_ALLELE);
                    for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
                        if (AlignmentUtils.isPartiallyEqual(a.getBase(taxon, site),badGeno)) {
                            if (a.isHeterozygous(taxon, site)==false) mna.setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                            else {
                                byte[] all= a.getBaseArray(taxon, site);
                                all[0]= (all[0]==NucleotideAlignmentConstants.INSERT_ALLELE||all[0]==NucleotideAlignmentConstants.GAP_ALLELE)?
                                        Alignment.UNKNOWN_ALLELE:Alignment.UNKNOWN_ALLELE;
                                mna.setBase(taxon, site, AlignmentUtils.getDiploidValue(all[0], all[1]));
                                }
                        }
                    }
                }
            }
            FilterAlignment fa= FilterAlignment.getInstance((a.retainsRareAlleles()==true)?mna:a, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));
            ExportUtils.writeToVCF(fa, file.substring(0, file.indexOf(inFiles[0].substring(inFiles[0].length()-6)))+"NoIndels.vcf.gz", '\t');
        }
    }
    
    public static void removeIndelsForBeagle(String dir, String fileType) {
        File[] inFiles= new File(dir).listFiles();
        
        for (File file:inFiles) {
            String readFile= null;
            if (file.isFile() && file.getName().endsWith(fileType)) readFile= file.getAbsolutePath();
            if (readFile==null) continue;
            Alignment a= ImportUtils.readGuessFormat(readFile, true);
            MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
            ArrayList<Integer> keepSites= new ArrayList<>();
            for (int site = 0; site < a.getSiteCount(); site++) {
                if (a.getMajorAllele(site)== NucleotideAlignmentConstants.GAP_ALLELE||
                        a.getMajorAllele(site)== NucleotideAlignmentConstants.INSERT_ALLELE||
                        a.getMinorAllele(site)== NucleotideAlignmentConstants.GAP_ALLELE||
                        a.getMinorAllele(site)== NucleotideAlignmentConstants.INSERT_ALLELE)
                    continue;
                keepSites.add(site);
                if (a.getAlleles(site).length>2) {
                    byte badGeno= AlignmentUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, NucleotideAlignmentConstants.INSERT_ALLELE);
                    for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
                        if (AlignmentUtils.isPartiallyEqual(a.getBase(taxon, site),badGeno)) {
                            mna.setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                        }
                    }
                }
            }
            FilterAlignment fa= FilterAlignment.getInstance(mna, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));
            ExportUtils.writeToVCF(fa, readFile.substring(0, readFile.indexOf(fileType))+"NoIndels.vcf.gz", '\t');
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
        
        dir= "/home/kls283/Documents/Imputation/beagle/";
        String[] files= new String[] {dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy282Allchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy282Allchr8.vcf.gz"};
        removeIndelsForBeagle(dir, ".hmp.h5");
        removeIndelsForBeagle(files);
        
        }
}