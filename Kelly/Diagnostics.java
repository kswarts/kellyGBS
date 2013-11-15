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
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
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
            outStream.close();
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
    
    public static void removeIndelsForBeagle(String dir, String fileType, boolean onlyPoly, double minSiteCov, double minTaxaCov) {
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
                if ((a.getTotalNotMissing(site)/a.getSequenceCount())<minSiteCov) continue;
                if (onlyPoly&&a.getMinorAlleleFrequency(site)<.0000000000000000000000001) continue;
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
            FilterAlignment fa= FilterAlignment.getInstance(mna, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));//filter sites
            ArrayList<Identifier> ids= new ArrayList<>();
            for (int taxon = 0; taxon < fa.getSequenceCount(); taxon++) {
                if (fa.getTotalNotMissingForTaxon(taxon)/fa.getSiteCount()<minTaxaCov) continue;
                ids.add(fa.getIdGroup().getIdentifier(taxon));
            }
            Alignment fat= FilterAlignment.getInstance(fa, new SimpleIdGroup(ids));//filter taxa
            ExportUtils.writeToVCF(fat, readFile.substring(0, readFile.indexOf(fileType))+"NoIndelsMinTCov"+minTaxaCov+"MinSCov"+minSiteCov+((onlyPoly==true)?"Poly":"Mono")+".vcf.gz", '\t');
            ExportUtils.writeToHDF5(fat, readFile.substring(0, readFile.indexOf(fileType))+"NoIndelsMinTCov"+minTaxaCov+"MinSCov"+minSiteCov+((onlyPoly==true)?"Poly":"Mono")+".hmp.h5");
        }
    }
    
    //match sites should be a subset of mod sites
    //this was written to subset the large file to make donors that match the specific populations to impute
    public static void matchSites(String dirToMatch, String fileType, String fileToMod) {
        Alignment mod= ImportUtils.readGuessFormat(fileToMod, true);
        File[] inFiles= new File(dirToMatch).listFiles();
        
        for (File file:inFiles) {
            String readFile= null;
            if (file.isFile() && file.getName().endsWith(fileType)) readFile= file.getAbsolutePath();
            if (readFile==null) continue;
            Alignment match= ImportUtils.readGuessFormat(readFile, true);
            ArrayList<Integer> keepModSites= new ArrayList<>();//keep sites in mod that have the same physical position in match
            ArrayList<Integer> keepMatchSites= new ArrayList<>();//remove sites in match that don't exist in mod
            int modSite= 0;
            for (int site = 0; site < match.getSiteCount(); site++) {
                modSite= mod.getSiteOfPhysicalPosition(match.getPositionInLocus(site), mod.getLocus(match.getLocusName(site)));
                if (modSite>-1) {keepModSites.add(modSite); keepMatchSites.add(site);}
            }
            FilterAlignment faMod= FilterAlignment.getInstance(mod, ArrayUtils.toPrimitive(keepModSites.toArray(new Integer[keepModSites.size()])));//filter sites
            ExportUtils.writeToHDF5(faMod, fileToMod.substring(0, fileToMod.length()-7)+"Match"+file.getName().substring(0, file.getName().length()-7)+".hmp.h5");
            if (keepMatchSites.size()!=match.getSiteCount()) {
                FilterAlignment faMatch= FilterAlignment.getInstance(match, ArrayUtils.toPrimitive(keepMatchSites.toArray(new Integer[keepModSites.size()])));
                ExportUtils.writeToHDF5(faMatch, readFile.substring(0, readFile.indexOf(fileType))+"Match"+fileToMod.substring(fileToMod.lastIndexOf('/'), file.getName().length()-7)+".hmp.h5");
            }
        }
    }
    
    public static void getCovByTaxon(String file) {
            Alignment a= ImportUtils.readGuessFormat(file, false);
            try {
            File outputFile = new File(file.substring(0, file.indexOf(".hmp")) + "TaxonCovDiagnostic.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
                int notMiss= 0;
                for (int site = 0; site < a.getSiteCount(); site++) {
                    if (a.getBase(taxon, site)==Alignment.UNKNOWN_DIPLOID_ALLELE) notMiss++;
                }
                double cov= notMiss/(double)a.getSiteCount();
                if (taxon!=0) outStream.writeBytes("\n"+Double.toString(cov));
                else outStream.writeBytes(Double.toString(cov));
            }
            outStream.close();
        }
        catch (Exception e) {
            System.out.println("Problem outputting file");
        }
    }
    
    public static void main(String[] args) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        String dir;
        
        dir= "/home/kls283/Documents/Imputation/";
//        String dir= "/Users/kls283/Desktop/Imputation/";
        String keyFile= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7.hmp.h5";
        String unmasked= dir+"AllZeaGBS_v2.7wDepth.hmp.h5";
//        maskedSitesDepth(keyFile, unmasked);
        
        dir= "/home/kls283/Documents/Imputation/beagle/";
        String[] files= new String[] {dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth5_Denom11StrictSubsetBy282Allchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy12S_RIMMA_Spanchr8.vcf.gz",
        dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy282Allchr8.vcf.gz"};
//        removeIndelsForBeagle(dir, ".hmp.h5", true, .1,.1);
//        removeIndelsForBeagle(files);
        
        dir= "/Users/kls283/Desktop/Imputation/";
        String[] covfiles= new String[] {dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy282Allchr8.hmp.h5",
            dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy12S_RIMMA_Spanchr8.hmp.h5",
            dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy282Allchr8.hmp.h5"};
//        for (String file:covfiles) {getCovByTaxon(file);}
        
        String file= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetByAmesTemperatechr8.hmp.h5";
        getCovByTaxon(file);
        
        }
}