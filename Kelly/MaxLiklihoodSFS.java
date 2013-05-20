/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import net.maizegenetics.pal.alignment.Alignment;

/**
 * Implements and maximum likelihood approach for reconstructing the 
 * site frequency spectrum proposed in Keightley and Halligan (Genetics
 * 2011) Takes a vcf with read depth information, fixed sites,>2 alleles,
 * and a designated outgroup.
 * @author kls283
 */
public class MaxLiklihoodSFS {
    public static String dir;
    public static byte diploidN= (byte) 0xff;
    
    public static int[] assignReadCodes(Alignment a, int refIndex, int siteIndex) {
        int[] codes= new int[a.getSequenceCount()];
       
        return codes;
    }
}
