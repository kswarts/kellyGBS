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
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author kls283
 */
public class PhaseHets {
    //Run Ed's code to get haplotypes for whole dataset
    //Only keep those with no hets
    //Use haplotypes to phase taxa (store two strands in different alignments)
    //Markov chain given nearest site as prior state. start with first present allele on one alignment, add alleles to alignment or alternate alignment based on conditional probability of seeing the states together,based on haplotypes
    public static final byte N= (byte) 0xFF;
    public static final byte A = (byte) 0x00;
    public static final byte C = (byte) 0x11;
    public static final byte G = (byte) 0x22;
    public static final byte T = (byte) 0x33;
    public static final byte INSERT = (byte) 0x44;
    public static final byte GAP = (byte) 0x55;
    
    public static void phase(String hapFile, String fileToPhase, String outFile) {
        Alignment haps= ImportUtils.readFromHapmap(hapFile, null);
        Alignment phase= ImportUtils.readFromHapmap(fileToPhase,null);
        if (haps.getSiteCount()!=phase.getSiteCount()) System.out.println("Sites do not match between alignment and haplotypes");
        Alignment[] doubleAlign= {phase,phase};
        MutableSingleEncodeAlignment gametes= MutableSingleEncodeAlignment.getInstance(doubleAlign);
        int plus= phase.getSequenceCount();
        int lastSite= 0;
        for (int taxon= 0; taxon<phase.getSequenceCount();taxon++) {
            for (int site= 0; site<phase.getSiteCount(); site++) {
                if (phase.getBase(taxon, site)!=N) {
                    if (lastSite==0) {//for the first not missing site the algorithm comes across
                        if (phase.isHeterozygous(taxon, site)==true) {//if het fill in appropriate diploid bases for both versions of the alignment
                            gametes.setBase(taxon, site, PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]));
                            gametes.setBase(taxon+plus, site, PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]));
                        }
                        else gametes.setBase(taxon+plus, site, N);//set second alignment to missing
                        lastSite= site;
                    }
                    else {
                        if (gametes.getBase(taxon+plus, lastSite)==N) {//if only working on first alignment
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                if (Math.random() <= getConditional(haps,site,firstBase,lastSite,phase.getBase(taxon, lastSite))) {
                                    gametes.setBase(taxon, site, firstBase);
                                    gametes.setBase(taxon+plus, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon, site, otherBase);
                                    gametes.setBase(taxon+plus, site, firstBase);
                                }
                            }
                            else {
                                if (Math.random() <= getConditional(haps,site,phase.getBase(taxon, site),lastSite,phase.getBase(taxon, lastSite))) gametes.setBase(taxon+plus, site, N);
                                else gametes.setBase(taxon, site, N);
                            }
                        }
                        else if (gametes.getBase(taxon, lastSite)==N) {//if only working on second alignment
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                if (Math.random() <= getConditional(haps,site,firstBase,lastSite,phase.getBase(taxon, lastSite))) {
                                    gametes.setBase(taxon+plus, site, firstBase);
                                    gametes.setBase(taxon, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, otherBase);
                                    gametes.setBase(taxon, site, firstBase);
                                }
                            }
                            else {
                                if (Math.random() <= getConditional(haps,site,phase.getBase(taxon, site),lastSite,phase.getBase(taxon, lastSite))) gametes.setBase(taxon, site, N);
                                else gametes.setBase(taxon+plus, site, N);
                            }
                        }
                        else {//if both alignments have alleles present at the previous site
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                double probFirstAlignOne= getConditional(haps,site,firstBase,lastSite,gametes.getBase(taxon, lastSite));
                                double probFirstAlignTwo= getConditional(haps,site,firstBase,lastSite,gametes.getBase(taxon+plus, lastSite));
                                if (probFirstAlignOne>=probFirstAlignTwo) {//assign the first base to the alignment with the highest cond probability
                                    gametes.setBase(taxon, site, firstBase);
                                    gametes.setBase(taxon+plus, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, firstBase);
                                    gametes.setBase(taxon, site, otherBase);
                                }
                            }
                            else {
                                byte currBase= phase.getBase(taxon, site);
                                double probAlignOne= getConditional(haps,site,currBase,lastSite,gametes.getBase(taxon, lastSite));
                                double probAlignTwo= getConditional(haps,site,currBase,lastSite,gametes.getBase(taxon+plus, lastSite));
                                if (probAlignOne>=probAlignTwo) {//assign the base to the alignment with the highest conditional probability and set the other alignment to missing at this site
                                    gametes.setBase(taxon, site, currBase);
                                    gametes.setBase(taxon+plus, site, N);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, currBase);
                                    gametes.setBase(taxon, site, N);
                                }
                            }
                        }
                        lastSite= site;
                    }
                }
            }
        }
        gametes.clean();
        ExportUtils.writeToHapmap(gametes, true, outFile, '\t', null);
    }
    
    public static double getConditional(Alignment haps, int focalSite, byte focalBase, int conditionalSite, byte conditionalBase) {
        double numCond= 0.0;
        double match= 0.0;
        for (int hap= 0;hap<haps.getSequenceCount();hap++) {
            if (haps.getBase(hap, conditionalSite)==conditionalBase) {
                numCond++;
                if (haps.getBase(hap, focalSite)==focalBase) match++;
            }
        }
        double cond= match/numCond!=0?match/numCond:0;
        
        return cond;
    }
    
    public static byte getDiploidBase(byte haploidBase) {
        byte diploidBase;
        if (haploidBase==(byte) 0xF) diploidBase= N;
        else if (haploidBase==(byte) 0x0) diploidBase= A;
        else if (haploidBase==(byte) 0x1) diploidBase= C;
        else if (haploidBase==(byte) 0x2) diploidBase= G;
        else if (haploidBase==(byte) 0x3) diploidBase= T;
        else if (haploidBase==(byte) 0x4) diploidBase= GAP;
        else diploidBase= INSERT;
        return diploidBase;
    }
}
