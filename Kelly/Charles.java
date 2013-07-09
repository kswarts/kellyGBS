/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;

/**
 *
 * @author kelly
 */
public class Charles {
    
    public static void myMethod(String mnah5FileName, String outSamplesFile, int someCutoff){
        MutableNucleotideAlignmentHDF5 mnah5 = MutableNucleotideAlignmentHDF5.getInstance(mnah5FileName);
        byte[] majAll= new byte[100];
        try{
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outSamplesFile), 655360));
          outStream.writeBytes("Major alleles for first 100 sites");
          for (int site = 0; site < 100; site++) {
              majAll[site]= mnah5.getMajorAllele(site);
              System.out.println(mnah5.getMajorAlleleAsString(site));
              outStream.writeBytes(mnah5.getMajorAlleleAsString(site));
          }
          outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
    }
    
    public static void main(String[] args) {
        String dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
        String hdf5= dir+"RIMMA_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1.hmp.txt.gz";
        myMethod(hdf5,"outFile",10);
        
    }
    
}