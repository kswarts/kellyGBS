/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author kelly
 * Imputes a file randomly based on MAF and proportion heterozygosity
 */
public class BaseImputation {
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    public static String dir;
    public static String baseFileName;
    
    public static void randomImpute(Alignment masked) {
        double[] hetFreq= new double[masked.getSequenceCount()];
        for (int taxon= 0;taxon<masked.getSequenceCount();taxon++) {
            hetFreq[taxon]= (double) masked.getHeterozygousCountForTaxon(taxon)/(double) masked.getTotalNotMissingForTaxon(taxon);
        }
        MutableSingleEncodeAlignment impute= MutableSingleEncodeAlignment.getInstance(masked, masked.getSequenceCount(), masked.getSiteCount());
        for (int site= 0;site<masked.getSiteCount();site++) {
            double minFreq= masked.getMinorAlleleFrequency(site);
            byte maj= masked.getMajorAllele(site);
            byte min= masked.getMinorAllele(site);
            String hetString= masked.getBaseAsString(site, maj)+masked.getBaseAsString(site, min);
            byte het= NucleotideAlignmentConstants.getNucleotideDiploidByte(hetString);
            for (int taxon= 0;taxon<masked.getSequenceCount();taxon++) {
                if (masked.getBase(taxon, site)==diploidN){
                    if (Math.random() <= minFreq) impute.setBase(taxon, site, min);
                    else impute.setBase(taxon, site, maj);
                    if (Math.random() <= hetFreq[taxon]) impute.setBase(taxon, site, het);
                }
            }
        }
        impute.clean();
        ExportUtils.writeToHapmap(impute, true, dir+baseFileName+"_imputed.hmp.txt", '\t', null);
    }
    
    public static void main (String[] args) {
        dir= "/Users/kelly/Documents/GBS/GroupImputation/";
        baseFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575";
        String maskedFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575_masked";
        BaseImputation.randomImpute(ImportUtils.readFromHapmap(dir+maskedFileName+".hmp.txt", true, null));
        
    }
}
