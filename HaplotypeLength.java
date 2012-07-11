/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;

/**Holds methods for comparing IBS and haplotype length between taxa and within groups of taxa
 *
 * @author kelly
 */
public class HaplotypeLength {
    public static String fileID;
    public static String chrNum;
    public static String[] chr= {"1","2","3","4","5","6","7","8","9","10"};

    /**records IBS between taxa and a reference (with each taxon as the reference, output to separate files) based on number of sites specified in the constructor**/
    public static void GetIBSBySite(int length) {

        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/NAM.hmp";
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,"1");
        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount());
        int remainder= inHapmap.getSiteCount()%length;
        System.out.println("remainder is: "+remainder);
        int lociByBlock= remainder==0?inHapmap.getSiteCount()/length:(inHapmap.getSiteCount()/length)+1;
        double[][][] IBSMatrix= new double[inHapmap.getSequenceCount()][inHapmap.getSequenceCount()][lociByBlock];
        double similarity= 0;
        int noScore= 0;
        char x;
        char y;

        for (int refTaxon= 0; refTaxon < inHapmap.getSequenceCount(); refTaxon++) {
            File newFile= new File("/Users/kelly/Documents/GBS/Ames/refTaxon_"+length+"_"+inHapmap.getFullTaxaName(refTaxon));

            try {
            newFile.createNewFile();
            newFile.setWritable(true);
            BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));
            
            for (int compareTaxon= 0; compareTaxon < inHapmap.getSequenceCount(); compareTaxon++) {
                BW.write(inHapmap.getTaxaName(compareTaxon)+"\t");
                
                for (int oneKBlock= 0; oneKBlock < lociByBlock-1; oneKBlock++) {//for the majority of the chromosome
                    for (int site= 0; site < length; site++) {
                        x= inHapmap.getBaseChar(refTaxon,(oneKBlock*length)+site);
                        y= inHapmap.getBaseChar(compareTaxon, (oneKBlock*length)+site);
                        if (x=='N'||y=='N') noScore+= 1;
                        else if (x=='R'||y=='R') {
                            if (x=='R'&&y=='R') similarity+=1.0;
                            else if (x=='R'&&(y=='A'||y=='G')) similarity+=0.5;
                            else if ((x=='A'||y=='G')&&y=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='Y'||y=='Y') {
                            if (x=='Y'&&y=='Y') similarity+=1.0;
                            else if(x == 'Y' && (y == 'C' || y == 'T')) similarity += 0.5;
                            else if ((x=='C'||y=='T')&&y=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='S'||y=='S') {
                            if (x=='S'&&y=='S') similarity+=1.0;
                            else if (x=='S'&&(y=='C'||y=='G')) similarity+=0.5;
                            else if ((x=='C'||y=='G')&&y=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='W'||y=='W') {
                            if (x=='W'&&y=='W') similarity+=1.0;
                            else if (x=='W'&&(y=='A'||y=='T')) similarity+=0.5;
                            else if ((x=='A'||y=='T')&&y=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='K'||y=='K') {
                            if (x=='K'&&y=='K') similarity+=1.0;
                            else if (x=='K'&&(y=='T'||y=='G')) similarity+=0.5;
                            else if ((x=='T'||y=='G')&&y=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='M'||y=='M') {
                            if (x=='M'&&y=='M') similarity+=1.0;
                            else if (x=='M'&&(y=='A'||y=='C')) similarity+=0.5;
                            else if ((x=='A'||y=='C')&&y=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='B'||y=='B') {
                            if (x=='B'&&y=='B') similarity+=1.0;
                            else if (x=='B'&&y!='A') similarity+=0.333;
                            else if (x!='A'&&y=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='D'||y=='D') {
                            if (x=='D'&&y=='D') similarity+=1.0;
                            else if (x=='D'&&y!='C') similarity+=0.333;
                            else if (x!='C'&&y=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='H'||y=='H') {
                            if (x=='H'&&y=='H') similarity+=1.0;
                            else if (x=='H'&&y!='G') similarity+=0.333;
                            else if (x!='G'&&y=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='V'||y=='V') {
                            if (x=='V'&&y=='V') similarity+=1.0;
                            else if (x=='V'&&y!='T') similarity+=0.333;
                            else if (x!='T'&&y=='V') similarity+=0.333;
                            else similarity+=0.0;
                        }
                        else if (x==y) similarity+=1.0;
                    }
                    IBSMatrix[refTaxon][compareTaxon][oneKBlock]= (similarity)/((double)length-(double)noScore);
                    System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t"+"Matrix: "+IBSMatrix[refTaxon][compareTaxon][oneKBlock]);
                    BW.write(Double.toString(IBSMatrix[refTaxon][compareTaxon][oneKBlock])+"\t");
                    similarity= 0.0;
                    noScore= 0;
                }

                if (lociByBlock!=0) {//for the remainder of the chromosome
                    for (int site= 0; site < remainder; site++) {
                        x= inHapmap.getBaseChar(refTaxon,((lociByBlock-1)*length)+site);
                        y= inHapmap.getBaseChar(compareTaxon, ((lociByBlock-1)*length)+site);
                        if (x=='N'||y=='N') noScore+= 1;
                        else if (x=='R'||y=='R') {
                            if (x=='R'&&y=='R') similarity+=1.0;
                            else if (x=='R'&&(y=='A'||y=='G')) similarity+=0.5;
                            else if ((x=='A'||y=='G')&&y=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='Y'||y=='Y') {
                            if (x=='Y'&&y=='Y') similarity+=1.0;
                            else if(x == 'Y' && (y == 'C' || y == 'T')) similarity += 0.5;
                            else if ((x=='C'||y=='T')&&y=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='S'||y=='S') {
                            if (x=='S'&&y=='S') similarity+=1.0;
                            else if (x=='S'&&(y=='C'||y=='G')) similarity+=0.5;
                            else if ((x=='C'||y=='G')&&y=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='W'||y=='W') {
                            if (x=='W'&&y=='W') similarity+=1.0;
                            else if (x=='W'&&(y=='A'||y=='T')) similarity+=0.5;
                            else if ((x=='A'||y=='T')&&y=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='K'||y=='K') {
                            if (x=='K'&&y=='K') similarity+=1.0;
                            else if (x=='K'&&(y=='T'||y=='G')) similarity+=0.5;
                            else if ((x=='T'||y=='G')&&y=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='M'||y=='M') {
                            if (x=='M'&&y=='M') similarity+=1.0;
                            else if (x=='M'&&(y=='A'||y=='C')) similarity+=0.5;
                            else if ((x=='A'||y=='C')&&y=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='B'||y=='B') {
                            if (x=='B'&&y=='B') similarity+=1.0;
                            else if (x=='B'&&y!='A') similarity+=0.333;
                            else if (x!='A'&&y=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='D'||y=='D') {
                            if (x=='D'&&y=='D') similarity+=1.0;
                            else if (x=='D'&&y!='C') similarity+=0.333;
                            else if (x!='C'&&y=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='H'||y=='H') {
                            if (x=='H'&&y=='H') similarity+=1.0;
                            else if (x=='H'&&y!='G') similarity+=0.333;
                            else if (x!='G'&&y=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='V'||y=='V') {
                            if (x=='V'&&y=='V') similarity+=1.0;
                            else if (x=='V'&&y!='T') similarity+=0.333;
                            else if (x!='T'&&y=='V') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x==y) similarity+=1.0;
                    }
                    IBSMatrix[refTaxon][compareTaxon][lociByBlock-1]= (similarity)/((double)remainder-(double)noScore);
//                    System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t"+"Matrix: "+IBSMatrix[refTaxon][compareTaxon][lociByThousandSites-1]);
                    BW.write(Double.toString(IBSMatrix[refTaxon][compareTaxon][lociByBlock-1])+"\t");
                    similarity= 0;
                    noScore= 0;
                }
                BW.write("\n");
            }
                BW.close();
                }
            catch (IOException e) {
                System.out.println(e);
            }
    }

}
    /**records proportion IBS for an entire chromosome**/
         public static void GetIBSForWholeChromosome() {

            String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/SS.hmp";
            Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,"1");
            System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount());
            double[][] IBSMatrix= new double[inHapmap.getSequenceCount()][inHapmap.getSequenceCount()];
            double similarity= 0.0;
            int noScore= 0;
            char x;
            char y;

            File newFile= new File("/Users/kelly/Documents/GBS/Ames/IBSWholeChromosome_Chr");

                try {
                newFile.createNewFile();
                newFile.setWritable(true);
                BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));
                BW.write("\t");
                for (int taxa= 0; taxa < inHapmap.getSequenceCount(); taxa++) {
                    BW.write(inHapmap.getTaxaName(taxa)+"\t");
                }
                BW.write("\n");

            for (int refTaxon= 0; refTaxon < inHapmap.getSequenceCount(); refTaxon++) {
                BW.write(inHapmap.getTaxaName(refTaxon)+"\t");

                for (int compareTaxon= 0; compareTaxon < inHapmap.getSequenceCount(); compareTaxon++) {
                                    
                    for (int site= 0; site < inHapmap.getSiteCount(); site++) {
                        x= inHapmap.getBaseChar(refTaxon,site);
                        y= inHapmap.getBaseChar(compareTaxon, site);
                        if (x=='N'||y=='N') noScore+= 1;
                        else if (x=='R'||y=='R') {
                            if (x=='R'&&y=='R') similarity+=1.0;
                            else if (x=='R'&&(y=='A'||y=='G')) similarity+=0.5;
                            else if ((x=='A'||y=='G')&&y=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='Y'||y=='Y') {
                            if (x=='Y'&&y=='Y') similarity+=1.0;
                            else if(x == 'Y' && (y == 'C' || y == 'T')) similarity += 0.5;
                            else if ((x=='C'||y=='T')&&y=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='S'||y=='S') {
                            if (x=='S'&&y=='S') similarity+=1.0;
                            else if (x=='S'&&(y=='C'||y=='G')) similarity+=0.5;
                            else if ((x=='C'||y=='G')&&y=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='W'||y=='W') {
                            if (x=='W'&&y=='W') similarity+=1.0;
                            else if (x=='W'&&(y=='A'||y=='T')) similarity+=0.5;
                            else if ((x=='A'||y=='T')&&y=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='K'||y=='K') {
                            if (x=='K'&&y=='K') similarity+=1.0;
                            else if (x=='K'&&(y=='T'||y=='G')) similarity+=0.5;
                            else if ((x=='T'||y=='G')&&y=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='M'||y=='M') {
                            if (x=='M'&&y=='M') similarity+=1.0;
                            else if (x=='M'&&(y=='A'||y=='C')) similarity+=0.5;
                            else if ((x=='A'||y=='C')&&y=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (x=='B'||y=='B') {
                            if (x=='B'&&y=='B') similarity+=1.0;
                            else if (x=='B'&&y!='A') similarity+=0.333;
                            else if (x!='A'&&y=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='D'||y=='D') {
                            if (x=='D'&&y=='D') similarity+=1.0;
                            else if (x=='D'&&y!='C') similarity+=0.333;
                            else if (x!='C'&&y=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='H'||y=='H') {
                            if (x=='H'&&y=='H') similarity+=1.0;
                            else if (x=='H'&&y!='G') similarity+=0.333;
                            else if (x!='G'&&y=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x=='V'||y=='V') {
                            if (x=='V'&&y=='V') similarity+=1.0;
                            else if (x=='V'&&y!='T') similarity+=0.333;
                            else if (x!='T'&&y=='V') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (x==y) similarity+=1.0;
//                        System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t");
                        
                    }
                        IBSMatrix[refTaxon][compareTaxon]= similarity/((double)inHapmap.getSiteCount()-(double)noScore);
//                        System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t"+"Matrix: "+IBSMatrix[refTaxon][compareTaxon]);
                        BW.write(Double.toString(IBSMatrix[refTaxon][compareTaxon])+"\t");
                        similarity= 0;
                        noScore= 0;
                    
                }
                    BW.write("\n");
                }
                    BW.close();
             }
                    
                catch (IOException e) {
                    System.out.println(e);
                }
    }   
        
/**like GetIBSBySite, but works based on a set physical position**/
    public static void GetIBSByPosition(int length) {

        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/SS.hmp";
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,"1");
        int lastPos= inHapmap.getPositionInLocus(inHapmap.getSiteCount()-1);
        int firstPos= inHapmap.getPositionInLocus(0);
        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount()+"\t"+"last Position: "+lastPos+"\t"+"first Position: "+firstPos);
        int remainder= lastPos%length;
        int numBlocks= remainder==0?lastPos/length:(lastPos/length)+1;
        System.out.println("numblocks is: "+numBlocks);
        double[][][] IBSMatrix= new double[inHapmap.getSequenceCount()][inHapmap.getSequenceCount()][numBlocks];
        double similarity= 0.0;
        int noScore= 0;
        int holdPosIncrease= 0;
        int sitesInBlock=0;
        char refAllele;
        char compareAllele;

        for (int refTaxon= 0; refTaxon < inHapmap.getSequenceCount(); refTaxon++) {
            File newFile= new File("/Users/kelly/Documents/GBS/Ames/refTaxon_"+inHapmap.getFullTaxaName(refTaxon));

            try {
            newFile.createNewFile();
            newFile.setWritable(true);
            BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));

            for (int compareTaxon= 0; compareTaxon < inHapmap.getSequenceCount(); compareTaxon++) {
                BW.write(inHapmap.getTaxaName(compareTaxon)+"\t");

                for (int block= 1; block <= numBlocks; block++) {

                    while (inHapmap.getPositionInLocus(holdPosIncrease+sitesInBlock) < (block!=numBlocks?(block*length):(lastPos-1))) {
                        refAllele= inHapmap.getBaseChar(refTaxon,holdPosIncrease+sitesInBlock);
                        compareAllele= inHapmap.getBaseChar(compareTaxon, holdPosIncrease+sitesInBlock);
                        sitesInBlock++;
                        if (refAllele=='N'||compareAllele=='N') noScore+= 1;
                        else   refAllele= inHapmap.getBaseChar(refTaxon,holdPosIncrease+sitesInBlock);
                        if (refAllele=='R'||compareAllele=='R') {
                            if (refAllele=='R'&&compareAllele=='R') similarity+=1.0;
                            else if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='G')&&compareAllele=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='Y'||compareAllele=='Y') {
                            if (refAllele=='Y'&&compareAllele=='Y') similarity+=1.0;
                            else if(refAllele == 'Y' && (compareAllele == 'C' || compareAllele == 'T')) similarity += 0.5;
                            else if ((refAllele=='C'||compareAllele=='T')&&compareAllele=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='S'||compareAllele=='S') {
                            if (refAllele=='S'&&compareAllele=='S') similarity+=1.0;
                            else if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='C'||compareAllele=='G')&&compareAllele=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='W'||compareAllele=='W') {
                            if (refAllele=='W'&&compareAllele=='W') similarity+=1.0;
                            else if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='T')&&compareAllele=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='K'||compareAllele=='K') {
                            if (refAllele=='K'&&compareAllele=='K') similarity+=1.0;
                            else if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='T'||compareAllele=='G')&&compareAllele=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='M'||compareAllele=='M') {
                            if (refAllele=='M'&&compareAllele=='M') similarity+=1.0;
                            else if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='C')&&compareAllele=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='B'||compareAllele=='B') {
                            if (refAllele=='B'&&compareAllele=='B') similarity+=1.0;
                            else if (refAllele=='B'&&compareAllele!='A') similarity+=0.333;
                            else if (refAllele!='A'&&compareAllele=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='D'||compareAllele=='D') {
                            if (refAllele=='D'&&compareAllele=='D') similarity+=1.0;
                            else if (refAllele=='D'&&compareAllele!='C') similarity+=0.333;
                            else if (refAllele!='C'&&compareAllele=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='H'||compareAllele=='H') {
                            if (refAllele=='H'&&compareAllele=='H') similarity+=1.0;
                            else if (refAllele=='H'&&compareAllele!='G') similarity+=0.333;
                            else if (refAllele!='G'&&compareAllele=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='V'||compareAllele=='V') {
                            if (refAllele=='V'&&compareAllele=='V') similarity+=1.0;
                            else if (refAllele=='V'&&compareAllele!='T') similarity+=0.333;
                            else if (refAllele!='T'&&compareAllele=='V') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele==compareAllele) similarity+=1.0;
                    }

                    IBSMatrix[refTaxon][compareTaxon][block-1]= (similarity)/((double)sitesInBlock-(double)noScore);
//                    System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t"+"Matrix: "+IBSMatrix[refTaxon][compareTaxon][block-1]);
                    BW.write(Double.toString(IBSMatrix[refTaxon][compareTaxon][block-1])+"\t");
                    similarity= 0;
                    noScore= 0;
                    holdPosIncrease+= sitesInBlock;
                    sitesInBlock= 0;                    
                }
                holdPosIncrease= 0;
                BW.write("\n");
            }
                BW.close();
                }
            catch (IOException e) {
                System.out.println(e);
            }
    }
}
    
    /**like GetIBSBySite, but works based on a set physical position. This one prints out one long output. Is faster than OneOutput2 because saves to memory and prints out at the end. Crashes with large files though because of memory issues**/
    public static void GetIBSByPositionOneOutput(int length) {

        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/"+fileID+".hmp";
        File newFile= new File("/Users/kelly/Documents/GBS/Ames/IBSByPos_"+fileID+"_chr"+chrNum+"_"+length+"_bp.txt");
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,chrNum);
        int lastPos= inHapmap.getPositionInLocus(inHapmap.getSiteCount()-1);
        int firstPos= inHapmap.getPositionInLocus(0);
        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount()+"\t"+"last Position: "+lastPos+"\t"+"first Position: "+firstPos);
        int remainder= lastPos%length;
        int numBlocks= remainder==0?lastPos/length:(lastPos/length)+1;
        System.out.println("numblocks is: "+numBlocks);
        double[][][] IBSMatrix= new double[inHapmap.getSequenceCount()][inHapmap.getSequenceCount()][numBlocks];
        double[][] IBSAvg= new double[numBlocks][4];
        double similarity= 0.0;
        int noScore= 0;
        int holdPosIncrease= 0;
        int sitesInBlock= 0;
        char refAllele;
        char compareAllele;

        for (int refTaxon= 0; refTaxon < inHapmap.getSequenceCount(); refTaxon++) {            
            for (int compareTaxon= refTaxon+1; compareTaxon < inHapmap.getSequenceCount(); compareTaxon++) {
                for (int block= 1; block <= numBlocks; block++) {
                    while (inHapmap.getPositionInLocus(holdPosIncrease+sitesInBlock) < (block!=numBlocks?(block*length):(lastPos))) {

                        refAllele= inHapmap.getBaseChar(refTaxon,holdPosIncrease+sitesInBlock);
                        compareAllele= inHapmap.getBaseChar(compareTaxon, holdPosIncrease+sitesInBlock);
                        sitesInBlock++;
                        if (refAllele=='N'||compareAllele=='N') noScore++;
                        else if (refAllele=='R'||compareAllele=='R') {
                            if (refAllele=='R'&&compareAllele=='R') similarity+=1.0;
                            else if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='G')&&compareAllele=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='Y'||compareAllele=='Y') {
                            if (refAllele=='Y'&&compareAllele=='Y') similarity+=1.0;
                            else if(refAllele == 'Y' && (compareAllele == 'C' || compareAllele == 'T')) similarity += 0.5;
                            else if ((refAllele=='C'||compareAllele=='T')&&compareAllele=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='S'||compareAllele=='S') {
                            if (refAllele=='S'&&compareAllele=='S') similarity+=1.0;
                            else if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='C'||compareAllele=='G')&&compareAllele=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='W'||compareAllele=='W') {
                            if (refAllele=='W'&&compareAllele=='W') similarity+=1.0;
                            else if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='T')&&compareAllele=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='K'||compareAllele=='K') {
                            if (refAllele=='K'&&compareAllele=='K') similarity+=1.0;
                            else if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='T'||compareAllele=='G')&&compareAllele=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='M'||compareAllele=='M') {
                            if (refAllele=='M'&&compareAllele=='M') similarity+=1.0;
                            else if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='C')&&compareAllele=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='B'||compareAllele=='B') {
                            if (refAllele=='B'&&compareAllele=='B') similarity+=1.0;
                            else if (refAllele=='B'&&compareAllele!='A') similarity+=0.333;
                            else if (refAllele!='A'&&compareAllele=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='D'||compareAllele=='D') {
                            if (refAllele=='D'&&compareAllele=='D') similarity+=1.0;
                            else if (refAllele=='D'&&compareAllele!='C') similarity+=0.333;
                            else if (refAllele!='C'&&compareAllele=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='H'||compareAllele=='H') {
                            if (refAllele=='H'&&compareAllele=='H') similarity+=1.0;
                            else if (refAllele=='H'&&compareAllele!='G') similarity+=0.333;
                            else if (refAllele!='G'&&compareAllele=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='V'||compareAllele=='V') {
                            if (refAllele=='V'&&compareAllele=='V') similarity+=1.0;
                            else if (refAllele=='V'&&compareAllele!='T') similarity+=0.333;
                            else if (refAllele!='T'&&compareAllele=='V') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele==compareAllele) similarity+=1.0;
                    }

                    IBSMatrix[refTaxon][compareTaxon][block-1]= (similarity)/((double)sitesInBlock-(double)noScore);
//                    System.out.println("similarity: "+similarity+"\t"+"noScore: "+noScore+"\t"+"Matrix: "+IBSMatrix[refTaxon][compareTaxon][block-1]);
                    IBSAvg[block-1][0]+= similarity;
                    IBSAvg[block-1][1]+= sitesInBlock;
                    IBSAvg[block-1][2]+= noScore;
                    holdPosIncrease+= sitesInBlock;
                    similarity= 0.0;
                    noScore= 0;                    
                    sitesInBlock= 0;                    
                }
                holdPosIncrease= 0;
                
            }
        }
        for (int block= 0; block < numBlocks; block++) {
            IBSAvg[block][3]= IBSAvg[block][0]/(IBSAvg[block][1]-IBSAvg[block][2]);
            System.out.println("Block "+block+":\t"+IBSAvg[block][3]+"\tSimilar:\t"+IBSAvg[block][0]+"\tSitesEval:\t"+IBSAvg[block][1]+"\tSitesNoScore:\t"+IBSAvg[block][2]);
        }

            try {
            newFile.createNewFile();
            newFile.setWritable(true);
            BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));
            
            BW.write("TaxonOne\tTaxonTwo");
            for (int g= 0; g<numBlocks; g++) {
                BW.write("\t"+(g+1));
            }
            BW.write("\n");
            for (int h=0; h < inHapmap.getSequenceCount(); h++) {
                for (int i=h+1; i < inHapmap.getSequenceCount(); i++) {
                    BW.write(inHapmap.getTaxaName(h)+"\t"+inHapmap.getTaxaName(i));
                    for(int j=0; j < numBlocks; j++) {
                        BW.write("\t"+Double.toString(IBSMatrix[h][i][j]));
                    } 
                BW.write("\n");
                }
            }
            BW.close();
            }
            catch (IOException e) {
                System.out.println(e);
            }
    
}

        /**like GetIBSBySite, but works based on a set physical position. Prints out one long output. Also includes sites evaluated  and missing data per block. Slower than OneOutput because prints to file as it goes - this saves memory space so should be able to use longer files**/
    public static void GetIBSByPositionOneOutput2(int length) {

        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/"+fileID+".hmp";
        File newFile= new File("/Users/kelly/Documents/GBS/Ames/IBSByPos2_"+fileID+"_chr"+chrNum+"_"+length+"_bp.txt");
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,chrNum);
        int lastPos= inHapmap.getPositionInLocus(inHapmap.getSiteCount()-1);
        int firstPos= inHapmap.getPositionInLocus(0);
        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount()+"\t"+"last Position: "+lastPos+"\t"+"first Position: "+firstPos);
        int remainder= lastPos%length;
        int numBlocks= remainder==0?lastPos/length:(lastPos/length)+1;
        System.out.println("numblocks is: "+numBlocks);
        double[] IBSBlockSim= new double[numBlocks];
        int[] IBSBlockNoScore= new int[numBlocks];
        int[] IBSBlockSites= new int[numBlocks];
        double[] IBSAvg= new double[numBlocks];
        double similarity= 0.0;
        int noScore= 0;
        int holdPosIncrease= 0;
        int sitesInBlock=0;
        char refAllele;
        char compareAllele;

        try {
            newFile.createNewFile();
            newFile.setWritable(true);
            BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));

            BW.write("TaxonOne\tTaxonTwo\tPosition\tSitesEvaluated\tMissingData\tPercentIBS\n");

        for (int refTaxon= 0; refTaxon < inHapmap.getSequenceCount(); refTaxon++) {
            for (int compareTaxon= 0; compareTaxon < inHapmap.getSequenceCount(); compareTaxon++) {
                for (int block= 1; block <= numBlocks; block++) {
                    while (inHapmap.getPositionInLocus(holdPosIncrease+sitesInBlock) < (block!=numBlocks?(block*length):(lastPos))) {

                        refAllele= inHapmap.getBaseChar(refTaxon,holdPosIncrease+sitesInBlock);
                        compareAllele= inHapmap.getBaseChar(compareTaxon, holdPosIncrease+sitesInBlock);
                        sitesInBlock++;
                        if (refAllele=='N'||compareAllele=='N') noScore+= 1;
                        else if (refAllele=='R'||compareAllele=='R') {
                            if (refAllele=='R'&&compareAllele=='R') similarity+=1.0;
                            else if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='G')&&compareAllele=='R') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='Y'||compareAllele=='Y') {
                            if (refAllele=='Y'&&compareAllele=='Y') similarity+=1.0;
                            else if(refAllele == 'Y' && (compareAllele == 'C' || compareAllele == 'T')) similarity += 0.5;
                            else if ((refAllele=='C'||compareAllele=='T')&&compareAllele=='Y') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='S'||compareAllele=='S') {
                            if (refAllele=='S'&&compareAllele=='S') similarity+=1.0;
                            else if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='C'||compareAllele=='G')&&compareAllele=='S') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='W'||compareAllele=='W') {
                            if (refAllele=='W'&&compareAllele=='W') similarity+=1.0;
                            else if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='T')&&compareAllele=='W') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='K'||compareAllele=='K') {
                            if (refAllele=='K'&&compareAllele=='K') similarity+=1.0;
                            else if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')) similarity+=0.5;
                            else if ((refAllele=='T'||compareAllele=='G')&&compareAllele=='K') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='M'||compareAllele=='M') {
                            if (refAllele=='M'&&compareAllele=='M') similarity+=1.0;
                            else if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')) similarity+=0.5;
                            else if ((refAllele=='A'||compareAllele=='C')&&compareAllele=='M') similarity+=0.5;
                            else similarity+=0;
                        }
                        else if (refAllele=='B'||compareAllele=='B') {
                            if (refAllele=='B'&&compareAllele=='B') similarity+=1.0;
                            else if (refAllele=='B'&&compareAllele!='A') similarity+=0.333;
                            else if (refAllele!='A'&&compareAllele=='K') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='D'||compareAllele=='D') {
                            if (refAllele=='D'&&compareAllele=='D') similarity+=1.0;
                            else if (refAllele=='D'&&compareAllele!='C') similarity+=0.333;
                            else if (refAllele!='C'&&compareAllele=='D') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='H'||compareAllele=='H') {
                            if (refAllele=='H'&&compareAllele=='H') similarity+=1.0;
                            else if (refAllele=='H'&&compareAllele!='G') similarity+=0.333;
                            else if (refAllele!='G'&&compareAllele=='H') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele=='V'||compareAllele=='V') {
                            if (refAllele=='V'&&compareAllele=='V') similarity+=1.0;
                            else if (refAllele=='V'&&compareAllele!='T') similarity+=0.333;
                            else if (refAllele!='T'&&compareAllele=='V') similarity+=0.333;
                            else similarity+=0;
                        }
                        else if (refAllele==compareAllele) similarity+=1.0;
                    }

                    if (refTaxon < compareTaxon) {
                    BW.write(inHapmap.getTaxaName(refTaxon)+"\t"+inHapmap.getTaxaName(compareTaxon)+"\t"+(block)+"\t"+(sitesInBlock)+"\t"+(noScore)+"\t"+((similarity)/((double)sitesInBlock-(double)noScore))+"\n");
                    }

                    IBSBlockSim[block-1]+= similarity;
                    IBSBlockNoScore[block-1]+= noScore;
                    IBSBlockSites[block-1]+= sitesInBlock;
                    similarity= 0;
                    noScore= 0;
                    holdPosIncrease+= sitesInBlock;
                    sitesInBlock= 0;
                }
                holdPosIncrease= 0;

            }
        }
        for (int a= 0; a < numBlocks; a++) {
            IBSAvg[a]= IBSBlockSim[a]/((double)IBSBlockSites[a]-(double)IBSBlockNoScore[a]);
            System.out.println("Block "+a+":\t"+IBSAvg[a]+"\tSimilar:\t"+IBSBlockSim[a]+"\tSitesEval:\t"+IBSBlockSites[a]+"\tSitesNoScore:\t"+IBSBlockNoScore[a]);
        }

            BW.close();
            }
            catch (IOException e) {
                System.out.println(e);
            }
}
    
    /**For characterizing haplotype diversity within a group. Compares
     * haplotype length between two randomly chosen taxa originating at a
     * randomly chosen site, iteratively. From the origin position taxa compared
     * based on site similarity to the right and left. Any site with missing
     * data in either taxa is skipped over (but recorded as missing).
     * Heterozygous SNP calls are considered similar if the corresponding
     * allele called is identically heterozygous or records either allele.
     * Haplotypes extend from one after the second diverged site to the left
     * (or the beginning of the chromosome) and extend to two before the
     * diverged site on the right (or the position of the last site). Skipping
     * over the first divergent site helps control for sequencing error. For
     * speed, works best with hapmap files filtered for non-polymorphic sites.
     * Output to tab-delimited text file.**/
    public static void GetAverageHaplotypeLength(int iterate) {
        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/"+fileID+".hmp";
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,chrNum);
        int lastSiteIndex= inHapmap.getSiteCount()-1;
        int lastPos= inHapmap.getPositionInLocus(lastSiteIndex);
        int firstPos= inHapmap.getPositionInLocus(0);
        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount()+"\t"+"last Position: "+lastPos+"\t"+"first Position: "+firstPos+"\t"+"times sampled: "+iterate);
        List<Integer> taxa= new ArrayList();
        List<Integer> sites= new ArrayList();
        int[][] hapLengths= new int[iterate][8]; //holds the taxa compared, the initial site, haplotype length, missing data, and "wobble" for each iteration
        long totalHapLength= 0;
        
        for (int i= 0; i < inHapmap.getSiteCount(); i++) { //initializes the arrayList that holds the indices for the sites
            sites.add(i);
        }
        
        for (int j= 0; j < inHapmap.getSequenceCount(); j++) { //initializes the arrayList that holds the indices for the taxa
            taxa.add(j);
        }
        
        for (int k= 0; k < iterate; k++) {
            
            Collections.shuffle(sites); //default uniform random shuffle of site indices
            Collections.shuffle(taxa); //default uniform random shuffle of taxa indices
            int taxonOne= taxa.get(0);
            int taxonTwo= taxa.get(1);
            int startPosIndex= sites.get(0);
            int currIncreasingPos= startPosIndex;
            int currDecreasingPos= startPosIndex-1; //minus one because starting position accounted for in first (increase) while loop
            int wobble=0; //records residual heterozygosity
            int sitesInHap= 0;
            int haplotype= 0; //holds length of haplotype by position
            int median= 0;
            int lowerBound= 0;
            int dissimilar= 0; //if alleles at a site differ
            int missingData= 0; //if either taxon has an 'N' at the site compared
            char refAllele;
            char compareAllele;            
            
            do { //compare alleles by site increasing from randomly chosen starting position between randomly chosen taxa
               refAllele= inHapmap.getBaseChar(taxonOne,currIncreasingPos);
               compareAllele= inHapmap.getBaseChar(taxonTwo, currIncreasingPos);
               
               if (refAllele=='N'||compareAllele=='N') {
                   missingData+= 1;
                   sitesInHap+= 1;
               }
               else if(refAllele == 'R' || compareAllele == 'R') {
                   if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')||compareAllele=='R') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='G')||refAllele=='R') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }

               else if (refAllele=='Y'||compareAllele=='Y') {
                   if (refAllele=='Y'&&(compareAllele=='C'||compareAllele=='T')||compareAllele=='Y') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='Y'&&(refAllele=='C'||refAllele=='T')||refAllele=='Y') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='S'||compareAllele=='S') {
                   if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')||compareAllele=='S') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='S'&&(refAllele=='C'||refAllele=='G')||refAllele=='S') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='W'||compareAllele=='W') {
                   if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')||compareAllele=='W') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='W'&&(refAllele=='A'||refAllele=='T')||refAllele=='W') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='K'||compareAllele=='K') {
                   if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')||compareAllele=='K') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='T'||refAllele=='G')||refAllele=='K') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='M'||compareAllele=='M') {
                   if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')||compareAllele=='M') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='C')||refAllele=='M') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }

               else if (refAllele=='B'||compareAllele=='B') {
                   if (refAllele=='B'&& compareAllele!='A') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='B'&& refAllele!='A') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }

               else if (refAllele=='D'||compareAllele=='D') {
                   if (refAllele=='D'&& compareAllele!='C') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='D'&& refAllele!='C') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }

               else if (refAllele=='H'||compareAllele=='H') {
                   if (refAllele=='H'&& compareAllele!='G') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='H'&& refAllele!='G') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='V'||compareAllele=='V') {
                   if (refAllele=='V'&& compareAllele!='T') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='V'&& refAllele!='T') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               else if (refAllele==compareAllele) sitesInHap+=1.0;
               else dissimilar += 1;

               currIncreasingPos++;
            } while (dissimilar < 2 && currIncreasingPos < lastSiteIndex);
            dissimilar= 0;
            
            do { //compare alleles by site decreasing from randomly chosen starting position between randomly chosen taxa
               refAllele= inHapmap.getBaseChar(taxonOne,currDecreasingPos);
               compareAllele= inHapmap.getBaseChar(taxonTwo, currDecreasingPos);
               
               if (refAllele=='N'||compareAllele=='N') {
                   missingData+= 1;
                   sitesInHap+= 1;
               }
               else if(refAllele == 'R' || compareAllele == 'R') {
                   if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')||compareAllele=='R') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='G')||refAllele=='R') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               
               else if (refAllele=='Y'||compareAllele=='Y') {
                   if (refAllele=='Y'&&(compareAllele=='C'||compareAllele=='T')||compareAllele=='Y') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='Y'&&(refAllele=='C'||refAllele=='T')||refAllele=='Y') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }
               
                else if (refAllele=='S'||compareAllele=='S') {
                   if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')||compareAllele=='S') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='S'&&(refAllele=='C'||refAllele=='G')||refAllele=='S') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }
               
                else if (refAllele=='W'||compareAllele=='W') {
                   if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')||compareAllele=='W') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='W'&&(refAllele=='A'||refAllele=='T')||refAllele=='W') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }

                else if (refAllele=='K'||compareAllele=='K') {
                   if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')||compareAllele=='K') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='T'||refAllele=='G')||refAllele=='K') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }
               
                else if (refAllele=='M'||compareAllele=='M') {
                   if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')||compareAllele=='M') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='C')||refAllele=='M') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               
               else if (refAllele=='B'||compareAllele=='B') {
                   if (refAllele=='B'&& compareAllele!='A') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='B'&& refAllele!='A') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               
               else if (refAllele=='D'||compareAllele=='D') {
                   if (refAllele=='D'&& compareAllele!='C') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='D'&& refAllele!='C') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               
               else if (refAllele=='H'||compareAllele=='H') {
                   if (refAllele=='H'&& compareAllele!='G') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='H'&& refAllele!='G') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
                }
               
                else if (refAllele=='V'||compareAllele=='V') {
                   if (refAllele=='V'&& compareAllele!='T') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else if (compareAllele=='V'&& refAllele!='T') {
                       sitesInHap+=1;
                       wobble+=1;
                   }
                   else dissimilar+=1;
               }
               else if (refAllele==compareAllele) sitesInHap+=1.0;
               else {
                   dissimilar += 1;
                   sitesInHap +=1;
               }
               
               currDecreasingPos--;
            } while (dissimilar < 2 && currDecreasingPos > 0);
            
            lowerBound= (inHapmap.getPositionInLocus(currDecreasingPos)+1);

            if (currDecreasingPos < 0 && currIncreasingPos > lastSiteIndex) {
                haplotype= lastPos;
                median= haplotype/2;
            }
            else if (currDecreasingPos < 0) {
                haplotype= (inHapmap.getPositionInLocus(currIncreasingPos)-1);//if the terminating polymorphic site is the first site on a chromosome, set to that
                median= haplotype/2;
            } 
            else if (currIncreasingPos > lastSiteIndex) {
                haplotype= lastPos-lowerBound;
                median= lowerBound+haplotype/2;
            }
            else {
                haplotype= (inHapmap.getPositionInLocus(currIncreasingPos)-1)-lowerBound;//minus one because these hapmaps only have polymorphic sites so length goes up to polymorphic site, covering space between markers
                median= lowerBound+haplotype/2;
            }
            
            hapLengths[k][0]= taxonOne;
            hapLengths[k][1]= taxonTwo;
            hapLengths[k][2]= inHapmap.getPositionInLocus(startPosIndex);
            hapLengths[k][3]= median;//the center of the haplotype
            hapLengths[k][4]= haplotype;
            hapLengths[k][5]= missingData;
            hapLengths[k][6]= wobble;
            hapLengths[k][7]= sitesInHap;
            totalHapLength+= haplotype;
        }
        System.out.println("Average Haplotype Length: "+"\t"+((double) totalHapLength/(double) iterate));

        File newFile= new File("/Users/kelly/Documents/GBS/Ames/avgHapLength"+fileID+"chr"+chrNum+iterate+"iterations");
        try {
                newFile.createNewFile();
                newFile.setWritable(true);
                BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));

                BW.write("taxonOne"+"\t"+"taxonTwo"+"\t"+"StartingPos"+"\t"+"Median"+"\t"+"HaplotypeLength"+"\t"+"SitesWithMissingData"+"\t"+"MatchingResidualHeterozygousSites"+"\t"+"SitesEvaluated"+"\n");
                for (int a= 0; a < iterate; a++) {
                    BW.write(inHapmap.getTaxaName(hapLengths[a][0])+"\t"+inHapmap.getTaxaName(hapLengths[a][1])+"\t"+hapLengths[a][2]+"\t"+hapLengths[a][3]+"\t"+hapLengths[a][4]+"\t"+hapLengths[a][5]+"\t"+hapLengths[a][6]+"\t"+hapLengths[a][7]+"\n");
                }
                
                BW.close();
        }
        catch (IOException e) {
                System.out.println(e);
            }
    }
    /**need same sites in reference hapmap and input hapmap**/
        public static void CompareHapToRef(String refHapmapName, String refTaxon) {
        String inHapmapFileName= "/Users/kelly/Documents/GBS/Ames/"+fileID+".hmp";
        String refHapmapFileName= "/Users/kelly/Documents/GBS/Ames/"+refHapmapName+".hmp";
        Alignment inHapmap= ImportUtils.readFromHapmap(inHapmapFileName,chrNum);
        Alignment refHapmap= ImportUtils.readFromHapmap(refHapmapFileName,chrNum);
        int lastSiteIndex= inHapmap.getSiteCount()-1;
        int lastPos= inHapmap.getPositionInLocus(lastSiteIndex);
        int firstPos= inHapmap.getPositionInLocus(0);
        int indexRefTaxon= refHapmap.getSequenceCount();
        int endPos= 0;
        int startPos= 0;
        int length= 0;

        for (int a=0; a < refHapmap.getSequenceCount(); a++) {
            if (refTaxon.equals(refHapmap.getTaxaName(a))) indexRefTaxon= a;
//            System.out.println(refHapmap.getTaxaName(a));
        }
        if (indexRefTaxon==refHapmap.getSequenceCount()) System.out.println("Need exact match for reference taxon name");

        System.out.println("No. of Taxa: "+inHapmap.getSequenceCount()+"\t"+"No of Sites: "+inHapmap.getSiteCount()+"\t"+"last Position: "+lastPos+"\t"+"first Position: "+firstPos+"\t"+"reference taxon: "+refTaxon); 
        char[] refStates=new char[refHapmap.getSiteCount()];
        int[][] output= new int[inHapmap.getSequenceCount()*inHapmap.getSiteCount()][9];//holds the taxa compared, the initial site, haplotype length, missing data, and "wobble" for each iteration
        int nextHap= 0;
        int wobble=0; //records residual heterozygosity
        int sitesInHap= 0;
        int dissimilar= 0; //if alleles at a site differ
        int missingData= 0; //if either taxon has an 'N' at the site compared

        for (int b= 0; b < refHapmap.getSiteCount(); b++) {
            refStates[b]= refHapmap.getBaseChar(indexRefTaxon, b);
        }
        if (refHapmap.getSiteCount()!=inHapmap.getSiteCount()) System.out.println("Need identical sites in reference and comparison hapmap files");


        for (int c= 0; c < inHapmap.getSequenceCount(); c++) {
            char refAllele;
            char compareAllele;         
            
            for (int k= 0;k < inHapmap.getSiteCount(); k++) {
            
                refAllele= refStates[k];
                compareAllele= inHapmap.getBaseChar(c, k);
                
                   if (refAllele=='N'||compareAllele=='N') {
                       missingData+= 1;
                       sitesInHap+= (sitesInHap<2)?0:1;
                   }

                   else if(refAllele == 'R' || compareAllele == 'R') {
                       if (refAllele=='R'&&(compareAllele=='A'||compareAllele=='G')||compareAllele=='R') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='G')||refAllele=='R') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                   }

                   else if (refAllele=='Y'||compareAllele=='Y') {
                       if (refAllele=='Y'&&(compareAllele=='C'||compareAllele=='T')||compareAllele=='Y') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='Y'&&(refAllele=='C'||refAllele=='T')||refAllele=='Y') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                    }

                    else if (refAllele=='S'||compareAllele=='S') {
                       if (refAllele=='S'&&(compareAllele=='C'||compareAllele=='G')||compareAllele=='S') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='S'&&(refAllele=='C'||refAllele=='G')||refAllele=='S') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                    }

                    else if (refAllele=='W'||compareAllele=='W') {
                       if (refAllele=='W'&&(compareAllele=='A'||compareAllele=='T')||compareAllele=='W') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='W'&&(refAllele=='A'||refAllele=='T')||refAllele=='W') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                    }

                    else if (refAllele=='K'||compareAllele=='K') {
                       if (refAllele=='K'&&(compareAllele=='T'||compareAllele=='G')||compareAllele=='K') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='R'&&(refAllele=='T'||refAllele=='G')||refAllele=='K') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                    }

                    else if (refAllele=='M'||compareAllele=='M') {
                       if (refAllele=='M'&&(compareAllele=='A'||compareAllele=='C')||compareAllele=='M') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='R'&&(refAllele=='A'||refAllele=='C')||refAllele=='M') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                   }

                   else if (refAllele=='B'||compareAllele=='B') {
                       if (refAllele=='B'&& compareAllele!='A') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='B'&& refAllele!='A') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                   }

                   else if (refAllele=='D'||compareAllele=='D') {
                       if (refAllele=='D'&& compareAllele!='C') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='D'&& refAllele!='C') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                   }

                   else if (refAllele=='H'||compareAllele=='H') {
                       if (refAllele=='H'&& compareAllele!='G') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='H'&& refAllele!='G') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                    }

                    else if (refAllele=='V'||compareAllele=='V') {
                       if (refAllele=='V'&& compareAllele!='T') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else if (compareAllele=='V'&& refAllele!='T') {
                           sitesInHap+=1;
                           wobble+=1;
                       }
                       else dissimilar+=1;
                   }
                   else if (refAllele==compareAllele) sitesInHap+=1.0;
                   else {
                    dissimilar += (sitesInHap<2)?0:1;
                    sitesInHap+= (sitesInHap<2)?0:1;
                   }

                if (dissimilar > 2||k==lastSiteIndex) {

                    startPos= inHapmap.getPositionInLocus(k-sitesInHap)+1;
                    endPos= inHapmap.getPositionInLocus(k)-1;
                    dissimilar= 0;
                    length= (startPos+endPos==0)?1:endPos-startPos;
                    
                    output[nextHap][0]= indexRefTaxon;
                    output[nextHap][1]= c;
                    output[nextHap][2]= startPos;
                    output[nextHap][3]= (startPos+length)/2;//the center of the haplotype
                    output[nextHap][4]= length;
                    output[nextHap][5]= missingData;
                    output[nextHap][6]= wobble;
                    output[nextHap][7]= sitesInHap;
                    output[nextHap][8]= nextHap+1;
                    
                    startPos= 0;
                    endPos= 0;
                    sitesInHap= 0;
                    wobble= 0;
                    missingData= 0;
                    nextHap++;
                }                
                }
            System.out.println("Average Haplotype Length Taxon"+inHapmap.getTaxaName(c)+":\t"+((double) lastPos/(double) nextHap));
                }
        

        File newFile= new File("/Users/kelly/Documents/GBS/Ames/compareHapLength2Mismatch"+fileID+"chr"+chrNum+"refTaxon_"+refTaxon);
        try {
                newFile.createNewFile();
                newFile.setWritable(true);
                BufferedWriter BW = new BufferedWriter(new FileWriter(newFile));

                BW.write("refTaxon"+"\t"+"compareTaxon"+"\t"+"startPos"+"\t"+"median"+"\t"+"hapLength"+"\t"+"missingData"+"\t"+"residHet"+"\t"+"sitesEvaluated"+"\n");
                int a= 0;
                do {
                    BW.write(refHapmap.getTaxaName(output[a][0])+"\t"+inHapmap.getTaxaName(output[a][1])+"\t"+output[a][2]+"\t"+output[a][3]+"\t"+output[a][4]+"\t"+output[a][5]+"\t"+output[a][6]+"\t"+output[a][7]+"\n");
                    a++;
                } while (output[a][8]!=0);

                BW.close();
        }
        catch (IOException e) {
                System.out.println(e);
            }
    }



    public static void main (String args[]) {
        fileID= "Ames_MAIZE_NoHets_602K_Chr10";
        chrNum= "10";
//        GetIBSBySite(100);
//        GetIBSByPosition(1000000);
        GetIBSByPositionOneOutput(1000000);
//        GetIBSByPositionOneOutput2(1000000);
//        GetIBSForWholeChromosome();
//        GetAverageHaplotypeLength(2000000);
//        CompareHapToRef("SS", "B73");

//        for(int a=0; a<10; a++) {
//            chrNum= chr[a];
//            GetAverageHaplotypeLength(40000);
//        }
//
        
    }
}
