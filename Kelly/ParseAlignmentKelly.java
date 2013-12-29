/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.*;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;
import org.apache.commons.lang.ArrayUtils;
import org.biojava3.core.util.Equals;

/**
 *
 * @author kelly
 */
public class ParseAlignmentKelly {
        public static String dir;
    public static byte diploidN= (byte) 0xff;
    
    //takes taxa from inFile and sites from refFile. Designed to compare 55k
   public static void checkSitesForIdentityByPosition(String inFile,String refFile, double siteThreshold, double taxonThreshold) {
       Alignment a= ImportUtils.readGuessFormat(dir+inFile, false);
       
       Alignment ref= ImportUtils.readGuessFormat(dir+refFile, true);
       MutableNucleotideAlignment mnaRef= MutableNucleotideAlignment.getInstance(ref);
       
       int[] refPos= ref.getPhysicalPositions();
       int[] matchTaxon= new int[a.getSequenceCount()];//holds the index of corresponding taxon in ref
       String[] refNames= new String[ref.getSequenceCount()];
       for (int taxon = 0; taxon < refNames.length; taxon++) {
           refNames[taxon]= ref.getIdGroup().getIdentifier(taxon).getNameLevel(0);
       }      
       for (int taxon = 0; taxon < matchTaxon.length; taxon++) {
           String unkName= a.getIdGroup().getIdentifier(taxon).getNameLevel(0);
           matchTaxon[taxon]= Arrays.binarySearch(refNames, unkName);
       }
       double[] badSites= new double[refPos.length];
       double[] halfSites= new double[refPos.length];
       double[] testedSites= new double[refPos.length];
       double[] badTaxa= new double[a.getSequenceCount()];
       double[] halfTaxa= new double[a.getSequenceCount()];
       double[] testedTaxa= new double[a.getSequenceCount()];
       for (int site= 0;site<a.getSiteCount();site++) {
           int refIndex= Arrays.binarySearch(refPos, a.getPositionInLocus(site));
           if (refIndex<0) continue;
           int refSite= ref.getSiteOfPhysicalPosition(a.getPositionInLocus(site), null);
           for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
               byte inBase= a.getBase(taxon, site);
               byte refBase= ref.getBase(matchTaxon[taxon], refSite);
               byte[] inBaseArray= a.getBaseArray(taxon, site);
               byte[] refBaseArray= ref.getBaseArray(matchTaxon[taxon], refSite);
               if (inBase==diploidN||refBase==diploidN) continue;
               testedSites[refIndex]+=1.0;
               testedTaxa[taxon]+=1.0;
               if (AlignmentUtils.isEqual(inBaseArray, refBaseArray)==true) continue;
               if (a.isHeterozygous(taxon, site)==true||ref.isHeterozygous(matchTaxon[taxon], refSite)==true) {
                   if (inBaseArray[0]==refBaseArray[0]||inBaseArray[1]==refBaseArray[0]||inBaseArray[0]==refBaseArray[1]||
                        inBaseArray[1]==refBaseArray[1]) {
                       halfSites[refIndex]+=1.0;
                       halfTaxa[taxon]+=1.0;
                   }
               }
               else {
                   badSites[refIndex]+=1.0;
                   badTaxa[taxon]+=1.0;
               }
           }
           if (badSites[refIndex]/testedSites[refIndex]>siteThreshold) mnaRef.clearSiteForRemoval(refIndex);
           System.out.println(a.getSNPID(site)+"\t"+badSites[refIndex]/testedSites[refIndex]
                   +"\t"+halfSites[refIndex]/testedSites[refIndex]);
       }
       boolean[] badIDs= new boolean[a.getSequenceCount()];
       for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
           System.out.println(a.getTaxaName(taxon)+"\t"+badTaxa[taxon]/testedTaxa[taxon]
                   +"\t"+halfTaxa[taxon]/testedTaxa[taxon]);
           if (badTaxa[taxon]/testedTaxa[taxon]>taxonThreshold) badIDs[taxon]= true;
           
       }
       IdGroup removeIDs= IdGroupUtils.idGroupSubset(a.getIdGroup(), badIDs);
       Alignment goodTaxa= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
       System.out.println("Write "+goodTaxa.getSequenceCount()+" taxa to file");
       ExportUtils.writeToHapmap(goodTaxa, true, dir+inFile.substring(0, inFile.lastIndexOf('h')-1) +"_cleanTaxa.hmp.txt.gz", '\t', null);
       
       mnaRef.clean();
       System.out.println("Write "+mnaRef.getSiteCount()+" sites to file");
       ExportUtils.writeToHapmap(mnaRef, true, dir+refFile.substring(0, refFile.lastIndexOf('h')-1)+"_cleanSites.hmp.txt.gz", '\t', null);
   }
   
   //this subsets an hdf5 based on taxa names taken from a text file (one taxon per line). Can be either permissive (where only the first
   //part of the name is compared) or strict, where only identical sample preps are selected. Print to hmp (h5==false) or hdf5 (h5==true)
   //if remove taxa equals true, the taxa listed in the text file will be removed and the rest returned, otherwise selected
   public static void subsetHDF5FromTxt(String dir, String inFileRoot, String taxaNamesRoot, boolean permissive, boolean removeTaxa, boolean keepDepth, boolean toVCF) {
       String[] nameArray= KellyUtils.readInTxtNames(dir+taxaNamesRoot+".txt", permissive);
       Alignment in= ImportUtils.readGuessFormat(dir+inFileRoot,true);
       if (!inFileRoot.endsWith(".h5")) {
           ExportUtils.writeToMutableHDF5(in, dir+inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+".hmp.h5");
           inFileRoot= inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+".hmp.h5";
       }
       MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot);
       mnah5.optimizeForTaxa(null);
       IdGroup h5IDs= mnah5.getIdGroup();
       if (permissive==true) System.out.println("permissive mode on (matches only the first part of the name before the first colon - will select duplicate samples of taxa)");
       if (removeTaxa==true) System.out.println("removeTaxa mode on. will remove taxa included in text file from input h5");
       if (keepDepth==true) System.out.println("output file will retain depth information");
       else System.out.println("output file will not retain depth information");
       ArrayList<Identifier> keep= new ArrayList<Identifier>();
       for (int h5Taxa = 0; h5Taxa < h5IDs.getIdCount(); h5Taxa++) {
           String currTaxon= permissive==true?h5IDs.getIdentifier(h5Taxa).getName():h5IDs.getIdentifier(h5Taxa).getFullName();
           if (Arrays.binarySearch(nameArray, currTaxon)>-1) keep.add(h5IDs.getIdentifier(h5Taxa));
       }
       SimpleIdGroup subIDs= new SimpleIdGroup(keep);
       Alignment subset;
       if (removeTaxa==false) {
           subset= FilterAlignment.getInstance(mnah5, subIDs);
           System.out.println("subsetting "+subset.getSequenceCount()+" taxa from "+h5IDs.getIdCount()+" taxa based on "+nameArray.length+" unique names");
           if (permissive==true) {ExportUtils.writeToMutableHDF5(subset, dir+inFileRoot+"PermissiveSubsetBy"+taxaNamesRoot+".hmp.h5", null, keepDepth);
           if (toVCF) {System.out.println("writing to vcf"); ExportUtils.writeToVCF(subset, dir+inFileRoot+"PermissiveSubsetBy"+taxaNamesRoot+".vcf.gz", '\t');}}
           else {ExportUtils.writeToMutableHDF5(subset, dir+inFileRoot+"StrictSubsetBy"+taxaNamesRoot+".hmp.h5", null, keepDepth);
           if (toVCF) {System.out.println("writing to vcf"); ExportUtils.writeToVCF(subset, dir+inFileRoot+"StrictSubsetBy"+taxaNamesRoot+".vcf.gz", '\t');}}
           System.out.println("writing to hdf5");
       }
       else {
           subset= FilterAlignment.getInstanceRemoveIDs(mnah5, subIDs);
           System.out.println("removing "+subIDs.getIdCount()+" taxa from "+mnah5.getSequenceCount()+" taxa based on "+nameArray.length+" unique names\n"+subset.getSequenceCount()+" taxa remain");
           if (permissive==true) {ExportUtils.writeToMutableHDF5(subset, dir+inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+"PermissiveExclusionBy"+taxaNamesRoot+".hmp.h5", null, keepDepth);
           if (toVCF) {System.out.println("writing to vcf"); ExportUtils.writeToVCF(subset, dir+inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+"PermissiveExclusionBy"+taxaNamesRoot+".vcf.gz", '\t');}}
           else {ExportUtils.writeToMutableHDF5(subset, dir+inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+"StrictExclusionBy"+taxaNamesRoot+".hmp.h5", null, keepDepth);
           System.out.println("writing to hdf5");
           if (toVCF) {System.out.println("writing to vcf"); ExportUtils.writeToVCF(subset, dir+inFileRoot.substring(0, inFileRoot.indexOf(".hmp"))+"PermissiveExclusionBy"+taxaNamesRoot+".vcf.gz", '\t');}}
           
       }
   }
   
   //will join alignments, keeping only those sites that are present in both
   public static void combineAlignments(String fileOne, String fileTwo) {
       Alignment a= ImportUtils.readGuessFormat(fileOne, true);
       Alignment b= ImportUtils.readGuessFormat(fileTwo, true);
       String outFile= fileOne.substring(0,fileOne.indexOf(".hmp"))+"_Join_"+fileTwo.substring(fileTwo.lastIndexOf('/'));
       int[] aPos= a.getPhysicalPositions();
       int[] bPos= b.getPhysicalPositions();
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a, a.getSequenceCount()+b.getSequenceCount(), aPos.length);
       boolean remove[]= new boolean[aPos.length];
       int newSize= 0;
       for (int site = 0; site < aPos.length; site++) {
           if (Arrays.binarySearch(bPos, aPos[site])>-1&&
                   (a.getMajorAllele(site)==b.getMajorAllele(site)||a.getMajorAllele(site)==b.getMinorAllele(site))&&
                   (b.getMajorAllele(site)==a.getMajorAllele(site)||b.getMajorAllele(site)==a.getMinorAllele(site))) {newSize++; continue;}
           else mna.clearSiteForRemoval(site);
       }
       IdGroup aTaxa= a.getIdGroup();
       for (int taxon = 0; taxon < b.getSequenceCount(); taxon++) {
           if (aTaxa.whichIdNumber(b.getIdGroup().getIdentifier(taxon))<0) continue;
           mna.addTaxon(b.getIdGroup().getIdentifier(taxon));
       }
       mna.clean();
       
       int[] newPos= mna.getPhysicalPositions();
       for (int taxon = 0; taxon < b.getSequenceCount(); taxon++) {
           int whichTaxon= mna.getIdGroup().whichIdNumber(b.getIdGroup().getIdentifier(taxon));
           if (whichTaxon<0) continue;
           byte[] newAlign= new byte[mna.getSiteCount()];
           for (int site = 0; site < b.getSiteCount(); site++) {
               int newSite= Arrays.binarySearch(newPos, bPos[site]);
               if (newSite < 0) continue;
               newAlign[newSite]= b.getBase(taxon, site);
           }
           mna.setBaseRange(whichTaxon, 0, newAlign);
       }
       mna.clean();
       ExportUtils.writeToHDF5(a, outFile);
   }
   
   public static void subsetAlignmentByFocusSites(String inFile, String outAdd, int[][] positions, int windowSize, boolean h5) {//the rows of positions should be physical position, the columns loci
       Alignment a= ImportUtils.readGuessFormat(inFile, true);
       ArrayList<Integer> all= new ArrayList<Integer>();
       for (int pos = 0; pos < positions.length; pos++) {
           ArrayList<Integer> sub= new ArrayList<Integer>();
           Locus currLocus= a.getLocus(String.valueOf(positions[pos][1]));
           int[] locusPos= Arrays.copyOfRange(a.getPhysicalPositions(), currLocus.getStart(),currLocus.getEnd());
           int closestSite= Math.abs(Arrays.binarySearch(locusPos, positions[pos][0]))+currLocus.getStart()-1;
           System.out.println("Closest site position to target site "+positions[pos][0]+" in locus "+positions[pos][1]+" is "+a.getPositionInLocus(closestSite)+" ("+closestSite+")");
           for (int site = Math.max(0,closestSite-windowSize); site < Math.min(closestSite+windowSize,a.getStartAndEndOfLocus(currLocus)[1]-1); site++) {
               sub.add(site);
               all.add(site);
           }
           if (h5==true) continue;
           else {
               int[] subSites=  ArrayUtils.toPrimitive(sub.toArray(new Integer[sub.size()]));
               Alignment subA= FilterAlignment.getInstance(a, subSites);
               if (subA.getLocus(0)!=currLocus) {
                   System.out.println("Locus not conveyed to filter alignment. Set new");
                   MutableNucleotideAlignment mnaSub= MutableNucleotideAlignment.getInstance(subA);
                   for (int site = 0; site < mnaSub.getSiteCount(); site++) {mnaSub.setLocusOfSite(site, currLocus);}
                   ExportUtils.writeToHapmap(mnaSub, true, inFile.substring(0, inFile.indexOf(".hmp"))+outAdd+"subString"+windowSize+"_c"+mnaSub.getLocusName(0)+".hmp.txt.gz", '\t', null);
                   continue;
               }
               ExportUtils.writeToHapmap(subA, true, inFile.substring(0, inFile.indexOf(".hmp"))+outAdd+"subString"+windowSize+"_c"+subA.getLocusName(0)+".hmp.txt.gz", '\t', null);
           }
       }
       if (h5==true) {
           int[] subSites=  ArrayUtils.toPrimitive(all.toArray(new Integer[all.size()]));
           ExportUtils.writeToMutableHDF5(a, inFile.substring(0, inFile.indexOf(".hmp"))+outAdd+"subString"+windowSize, subSites);
       }
   }
   
   //donor file root should contain ".gX." to represent the chromosome and segment ie c1s4
   public static void findSubsetDonor(String donorFile, String matchFile) {
       Alignment match= ImportUtils.readGuessFormat(matchFile, true);
       File theDF=new File(donorFile);
       String prefilter=theDF.getName().split(".gX.")[0]+".gc"; //grabs the left side of the file
       String prefilterOld=theDF.getName().split("s\\+")[0]+"s"; //grabs the left side of the file
       ArrayList<String> d=new ArrayList<String>();
       for (File file : theDF.getParentFile().listFiles()) {
            if(file.getName().equals(theDF.getName())) {d.add(file.toString());}
            if(file.getName().startsWith(prefilter)) {d.add(file.toString());}
            if(file.getName().startsWith(prefilterOld)) {d.add(file.toString());}
        }
        ArrayList<Alignment> aligns= new ArrayList<>();
        for (String tryFile:d){
            Alignment donorAlign=ImportUtils.readGuessFormat(tryFile, true);
            if (Arrays.binarySearch(match.getLoci(),donorAlign.getLoci()[0])<0) continue;
            else for (int pos:donorAlign.getPhysicalPositions()) {
                if (Arrays.binarySearch(match.getPhysicalPositions(),pos)<0) continue;
                else {
                    aligns.add(donorAlign);
                    break;
                }
            }
        }
        for (Alignment don:aligns) {
            String outFileName= donorFile.substring(0, donorFile.indexOf(".gX"))+"MatchTo"+matchFile.substring(matchFile.lastIndexOf('/'), matchFile.lastIndexOf(".hmp.txt"))+".gc"+don.getLoci()[0]+"s"+aligns.indexOf(don)+".hmp.txt";
            ArrayList<Integer> sub= new ArrayList<>();
            for (int matchPos:match.getPhysicalPositions()) {if (Arrays.binarySearch(don.getPhysicalPositions(), matchPos)>-1) sub.add(matchPos);}
            int[] subSites= ArrayUtils.toPrimitive(sub.toArray(new Integer[sub.size()]));
            Alignment outDonor= FilterAlignment.getInstance(don, subSites);
            if (Equals.equal(outDonor.getLoci(),don.getLoci())==false) {
                MutableNucleotideAlignment outMnaDonor= MutableNucleotideAlignment.getInstance(outDonor);
                for (int site = 0; site < outMnaDonor.getSiteCount(); site++) {outMnaDonor.setLocusOfSite(site, don.getLoci()[0]);}
                ExportUtils.writeToHapmap(outMnaDonor, true, outFileName, '\t', null);
                continue;
            }
            ExportUtils.writeToHapmap(outDonor, true, outFileName, '\t', null);
        }
        
   }
   //cuts up the genome into chunks and outputs delimited text with Major==1, Minor==2, Het==3. Missing==4
   public static void outputCodedTabTextForDong(String inFile, int chunkSize, int howManyInOne, String delimiter) {
       MutableNucleotideAlignmentHDF5 a= MutableNucleotideAlignmentHDF5.getInstance(inFile);
       int startSite= 0;
       for (int part = 1; part < Math.floor(a.getSequenceCount()/(chunkSize*howManyInOne))+1; part++) {
           int size= chunkSize*howManyInOne;
           try {
               BufferedWriter bw= new BufferedWriter(new FileWriter(inFile.substring(0, inFile.indexOf(".hmp.h5"))+"_part"+part+".txt"));
               for (int site = startSite; site < Math.min(part*size, a.getSiteCount()); site++) {bw.write(delimiter+site);}
               for (int taxa = 0; taxa < a.getSequenceCount(); taxa++) {
                   bw.write("\n"+a.getTaxaName(taxa));
                   for (int site = startSite; site < Math.min(part*size, a.getSiteCount()); site++) {
                       int siteCode= a.getBase(taxa, site)==diploidN?4:a.isHeterozygous(taxa, site)?3:a.getBaseArray(taxa, site)[0]==a.getMinorAllele(site)?2:1;
                       bw.write(delimiter+siteCode);
                   }
               }
               bw.flush();
               bw.close();
           }
           catch (Exception e) {
               System.out.println ("Problem with writing file "+part);
           }
           
           startSite+= size;
        }
   }
   
   public static void unphasedUnrelatedForBeagle(String inFileName) {
       Alignment a= ImportUtils.readGuessFormat(inFileName, true);
       Locus[] allLoci= a.getLoci();
       for (Locus currLocus:allLoci) {
           try {
                File outputFile= new File(inFileName.substring(0, inFileName.indexOf(".hmp"))+"Beagle_c"+currLocus.getName()+".txt");
                DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
                outStream.writeBytes("I\tid");
                for (int taxa = 0; taxa < a.getSequenceCount(); taxa++) {outStream.writeBytes("\t"+a.getFullTaxaName(taxa));}
                for (int site = currLocus.getStart(); site < currLocus.getEnd()+1; site++) {
                    outStream.writeBytes("\nM\tl"+a.getLocusName(site)+"p"+a.getPositionInLocus(site));
                    for (int taxa = 0; taxa < a.getSequenceCount(); taxa++) {
                        if (a.isHeterozygous(taxa, site)==true) outStream.writeBytes("\t"+a.getBaseAsStringArray(taxa, site)[0]+"/"+a.getBaseAsStringArray(taxa, site)[1]);
                        else outStream.writeBytes("\t"+a.getBaseAsString(taxa, site)+"/"+a.getBaseAsString(taxa, site));
                    }
                }
           }
           catch (Exception e) {
           System.out.println(e);
           }
       }
   }
   
   public static void extractPhaseFromImputedBeagle4 (String beagleFileName, String unimputedFileName) {
       //pull file statistics from beagle read in as vcf and set up new mna
       Alignment beagle= ImportUtils.readFromVCF(beagleFileName, null);
       Alignment unimputed= ImportUtils.readGuessFormat(unimputedFileName, true);
       if (Arrays.equals(beagle.getSNPIDs(),unimputed.getSNPIDs())==false) System.out.println("WARNING: Beagle and unimputed file to do not contain identical sites. Make sure to use an identical unimputed and imputed file");
       if (IdGroupUtils.isEqualIgnoringOrder(beagle.getIdGroup(),unimputed.getIdGroup())) System.out.println("WARNING: Beagle and unimputed file to do not contain the same taxa. Make sure to use an identical unimputed and imputed file");;
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(beagle, beagle.getSequenceCount()*2, beagle.getSiteCount());
       for (int taxon = 0; taxon < beagle.getSequenceCount(); taxon++) {//set new taxon names for phased haps
           String currID= beagle.getIdGroup().getIdentifier(taxon).getFullName();
           String IDOne= currID.substring(0, currID.indexOf(':'))+"Ref"+currID.substring(currID.indexOf(':'));
           String IDTwo= currID.substring(0, currID.indexOf(':'))+"Alt"+currID.substring(currID.indexOf(':'));
           mna.setTaxonName(taxon, new Identifier(IDOne));
           mna.setTaxonName(beagle.getSequenceCount()+taxon, new Identifier(IDTwo));
       }
       if (IdGroupUtils.isEqualIgnoringOrder(beagle.getIdGroup(), unimputed.getIdGroup())==false) System.out.println("Taxa names not equal. Ignore output, it's garbage");
       //go through and for each genotype that is not unknown in the unimputed, sort imputed result into ref/alt
       for (int site = 0; site < beagle.getSiteCount(); site++) {
           int unimpSite= unimputed.getSiteOfPhysicalPosition(beagle.getPositionInLocus(site), beagle.getLocus(site));
           for (int taxon = 0; taxon < beagle.getSequenceCount(); taxon++) {
               mna.setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);//wipe the current allele
               if (unimputed.getBase(taxon, unimpSite)==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
               if (AlignmentUtils.isPartiallyEqual(beagle.getBase(taxon, site), unimputed.getBase(taxon, unimpSite))==false) {
                   System.out.println("Do not partialy agree at position "+beagle.getPositionInLocus(site)); continue;
               }
               if (AlignmentUtils.isHeterozygous(beagle.getBase(taxon, site))==true) {
                   mna.setBase(taxon, site, AlignmentUtils.getDiploidValue(beagle.getBaseArray(taxon, site)[0], Alignment.UNKNOWN_ALLELE));
                   mna.setBase(taxon+beagle.getSequenceCount(), site, AlignmentUtils.getDiploidValue(beagle.getBaseArray(taxon, site)[1], Alignment.UNKNOWN_ALLELE));
               }
               else {
                   mna.setBase(taxon, site, AlignmentUtils.getDiploidValue(beagle.getBaseArray(taxon, site)[0], Alignment.UNKNOWN_ALLELE));
                   mna.setBase(taxon+beagle.getSequenceCount(), site, AlignmentUtils.getDiploidValue(beagle.getBaseArray(taxon, site)[0], Alignment.UNKNOWN_ALLELE));
               }
           }
       }
       mna.clean();
       ExportUtils.writeToHDF5(mna, unimputedFileName.substring(0, unimputedFileName.length()-7)+"BeagleHaps.hmp.h5");
   }
   
   public static void subsetChrKeepAllSites(String inFile, String chr, boolean vcf) {
       String outFile= inFile.substring(0, inFile.indexOf(".hmp"))+"chr"+chr+".hmp.h5";
       Alignment a= ImportUtils.readGuessFormat(inFile, true);
       int[] startEnd = a.getStartAndEndOfLocus(a.getLocus(chr));
       FilterAlignment fa= FilterAlignment.getInstance(a, startEnd[0], startEnd[1]-1);
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(fa);
       ExportUtils.writeToMutableHDF5(mna, outFile+".hmp.h5");
       if (vcf) ExportUtils.writeToVCF(mna, outFile+".vcf.gz", '\t');
   }
   
   public static void replaceGenotypes(String modFileName, String replaceDir, String replaceFileType) {
       TasselPrefs.putAlignmentRetainRareAlleles(false);
       File[] inFiles= new File(replaceDir).listFiles();
       ArrayList<Alignment> aligns= new ArrayList<>();
       for (File file:inFiles) {
           String readFile= null;
           if (file.isFile()&&file.getName().contains(replaceFileType)) readFile= file.getAbsolutePath();
           if (readFile==null) continue;
           aligns.add(ImportUtils.readGuessFormat(readFile, true));
       }
       if (aligns.isEmpty()) {System.out.println("No alignments to substitute in directory");return;}
       MutableNucleotideAlignmentHDF5 old= MutableNucleotideAlignmentHDF5.getInstance(modFileName);
       String newFile= modFileName.substring(0, modFileName.indexOf(".hmp.h5"))+"modFrom"+replaceDir.substring(replaceDir.lastIndexOf('/'))+".hmp.h5";
       ExportUtils.writeToMutableHDF5(old, newFile);
       MutableNucleotideAlignmentHDF5 mod= MutableNucleotideAlignmentHDF5.getInstance(newFile);
       Identifier[] modTaxa= new Identifier[mod.getSequenceCount()];
       for (int taxon = 0; taxon < modTaxa.length; taxon++) {modTaxa[taxon]= mod.getIdGroup().getIdentifier(taxon);}
       for(Alignment sub:aligns) {
           int[] currModRange= new int[]{mod.getSiteOfPhysicalPosition(sub.getPositionInLocus(0), sub.getLocus(0)),
               mod.getSiteOfPhysicalPosition(sub.getPositionInLocus(sub.getSiteCount()-1), sub.getLocus(sub.getSiteCount()-1))};//the start and end sites in the mod file for the curr alignment. the end site is inclusive
           TreeMap<Integer,Integer> whichSubSite= new TreeMap<>();
           int[] subSites= sub.getPhysicalPositions();
           for (int site = currModRange[0]; site < currModRange[1]+1; site++) {whichSubSite.put(site,Arrays.binarySearch(subSites, mod.getPositionInLocus(site)));}
           for (int taxon = 0; taxon < sub.getSequenceCount(); taxon++) {
               int modTaxon= Arrays.binarySearch(modTaxa, sub.getIdGroup().getIdentifier(taxon));
               if (modTaxon<0) {System.out.println("Taxon "+sub.getTaxaName(taxon)+" cannot be found in modFile");continue;}
               byte[] curr= mod.getBaseRow(modTaxon);
               for (int site = currModRange[0]; site < currModRange[1]+1; site++) {//site here is the site for the modFile
                   if (whichSubSite.get(site)<0) continue;
                   curr[site]= sub.getBase(taxon, whichSubSite.get(site));
               }
               mod.setAllBases(modTaxon, curr);
           }
       }
       mod.clean();
   }
   
   public static void main (String args[]) {
       TasselPrefs.putAlignmentRetainRareAlleles(false);
//       dir= "";
//       String inFileName= dir+"";
//       String refFileName= dir+"";
//       checkSitesForIdentityByPosition(inFileName,refFileName,.1,.1);
       
////       dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
//       dir= "/Users/kellyadm/Desktop/Imputation2.7/";
////       //String textToSubsetRoot= "Ames(no EP or GEM)";
////       String textToSubsetRoot= "12S_RIMMA_Span_SEED";
//       String textToSubsetRoot= "Ames380";
//       String inFileNameRoot= "AllZeaGBSv27";
////       //imputation subset for carotenoids (Ames inbreds and also 12S, RIMMA, spanish landraces)
//////       subsetHDF5FromTxt(inFileNameRoot,textToSubsetRoot,false,false, true);
////       //subset out landraces to replace with those inbred by HMM
////       subsetHDF5FromTxt(inFileNameRoot,textToSubsetRoot,true,true, false);
//       //subset 380 inbreds from Ames (randomized in excel)
//       subsetHDF5FromTxt(inFileNameRoot,textToSubsetRoot,false,false,false);
       
//       //filter alignment subsetting doesn't work for h5 mutable alignments
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
////       String[] inFileName= new String[] {dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._masked_Depth5_Denom11.hmp.h5",
////       dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._maskKey_Depth5_Denom11.hmp.h5",
////       dir+"AllZeaGBSv27StrictSubsetByAmes(no EP or GEM)._masked_Depth5_Denom11.hmp.h5",
////       dir+"AllZeaGBSv27StrictSubsetByAmes(no EP or GEM)._maskKey_Depth5_Denom11.hmp.h5"};
////       String[] inFileName=  new String[] {dir+"AllZeaGBS_v2.7InbredFor12S_RIMMA_Span_SEED.hmp.h5"};
//       String[] inFileName=  new String[] {dir+"AllZeaGBSv27.hmp.h5"};
//       int[][] carotenoid= new int[][]{{86833000,1},{44444500,2},{82019000,6},{138886000,8}};//Lut1,zeaxanthin epoxidase,Y1,lycE
//       for (String name:inFileName) {subsetAlignmentByFocusSites(name,"Carotenoid",carotenoid,2000,false);}
       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
//       String h5= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._masked_Depth5_Denom17.hmp.h5";
//       outputCodedTabTextForDong(h5, 3968, 10, "\t");
       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
//       String inFile= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._masked_Depth5_Denom11CarotenoidsubString1000.hmp.h5";
////       String inFile= "AllZeaGBSv27StrictSubsetByAmes(no EP or GEM)._masked_Depth5_Denom11CarotenoidsubString1000.hmp.h5";
//       unphasedUnrelatedForBeagle(inFile);
       
//       dir= "/Users/kellyadm/Desktop/Imputation2.7/";
//       String donorRoot= dir+"donors/AllZeaGBS_v2.7InbredFor12S_RIMMA_Span_SEED_HaplotypeMergeInbredLandrace8k";
//       String matchFile= dir+"carotenoid/AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._masked_Depth5_Denom11CarotenoidsubString2000_cX.hmp.txt.gz";
//       findSubsetDonor(donorRoot, matchFile);
       
//       dir= "/Users/kls283/Desktop/Imputation/";
////       dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
////       String[] inFile= new String[]{dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._masked_Depth5_Denom11.hmp.h5",
////           dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._maskKey_Depth5_Denom11.hmp.h5"};
//       String subset= "12S_RIMMA_Span";
//       subset= "282All";
////       subset= "Ames380";
//       String depth= "7";
//       String denom= "7";
//       String[] inFile= new String[]{dir+"AllZeaGBS_v2.7wDepth_masked_Depth"+depth+"_Denom"+denom+"StrictSubsetBy"+subset+".hmp.h5"};//,
//           //dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth"+depth+"_Denom"+denom+"StrictSubsetBy"+subset+".hmp.h5"};
//       String chr= "8";
//       for (String in:inFile) {subsetChrKeepAllSites(in, chr, true);}
       
       //subsets for runtime test
       dir= "/kls283/Documents/Imputation";
       String inFile= "AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetByLessThan.01Het.rand5kNoIndelsMinTCov0.1MinSCov0.1Poly.hmp.txt.gz";
       String taxa= "LessThan.01Het.rand3k"; 
//       subsetHDF5FromTxt(dir,inFile,taxa,false,false,false,true);
       
       //test replaceGenotypes
       dir= "/home/kls283/Documents/Imputation/";
       String modFile= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7.hmp.h5";
       String replaceDir= "/Users/kls283/Desktop/Imputation/beagle/new";
       String fileType= ".vcf.gz";
       replaceGenotypes(modFile,replaceDir,fileType);
   }
    
}
