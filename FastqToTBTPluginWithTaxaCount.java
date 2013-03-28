/*
 * fastqToTBTPlugin
 */


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.util.DirectoryCrawler;

// libraries needed for Plugin functionality
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import java.awt.Frame;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaShort;
import org.apache.log4j.Logger;

//libraries added by KLS
import java.util.Arrays;
import net.maizegenetics.util.ArgsEngine;

/**
 * This pipeline converts a series of fastq files to TagsByTaxa files (one per fastq file).
 * It requires a list of existing tags (Tags object), which may come from a TagCounts file or TOPM file.
 *
 * @author james
 */
public class FastqToTBTPluginWithTaxaCount extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FastqToTBTPluginWithTaxaCount.class);
    private ArgsEngine myArgsEngine = null;

    private String[] myFastqFileS = null;
    private String   myKeyFile = null;
    private String   myEnzyme = null;
    private String   myOutputDir = null;
    private int      myMinCount = 1;
    private Tags     myMasterTags = null;
    private static int maxGoodReads = 200000000; // maximum number of good barcoded reads expected in a fastq file
    private boolean useTBTByte = false;
    private boolean useTBTShort = false;

    public FastqToTBTPluginWithTaxaCount() {
        super(null, false);
    }

    public FastqToTBTPluginWithTaxaCount(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

//        if ((myEnzyme == null) || (myEnzyme.length() == 0)) {
//            printUsage();
//            throw new IllegalStateException("performFunction: enzyme must be set.");
//        }
        // TODO - More checks to validate parameters...

        matchTagsToTaxa(myFastqFileS, myKeyFile, myEnzyme, myMasterTags, myOutputDir, myMinCount, useTBTByte, useTBTShort);
        return null;
    }

    private void printUsage() {
        myLogger.info(
            "\nUsage is as follows:\n"
                + "-i  Input directory containing .fastq files\n"
                + "-k  Barcode key file\n"
                + "-e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
                + "-o  Output directory\n"
                + "-s  Max good reads per lane. (Optional. Default is 200,000,000).\n"
                + "-c  Minimum taxa count within a fastq file for a tag to be output (default 1)\n"  // Nb: using TagsByTaxaBit, so max count PER TAXON = 1
                + "-y  Output to tagsByTaxaByte (tag counts per taxon from 0 to 127) instead of tagsByTaxaBit (0 or 1)\n"
                + "One of either:\n"
                + "    -t  Tag count file, OR A\n"
                + "    -m  Physical map file containing alignments\n");
    }

    @Override
    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-k", "--key-file", true);
            myArgsEngine.add("-e", "--enzyme", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-s", "--max-reads", true);
            myArgsEngine.add("-c", "--min-count", true);
            myArgsEngine.add("-y", "--TBTByte", false);
            myArgsEngine.add("-sh", "--TBTShort", false);
            myArgsEngine.add("-t", "--tag-count", true);
            myArgsEngine.add("-m", "--physical-map", true);
        }
        myArgsEngine.parse(args);

        String tempDirectory = myArgsEngine.getString("-i");
        
        if(myArgsEngine.getBoolean("-s")){ maxGoodReads = Integer.parseInt(myArgsEngine.getString("-s"));}

        if (tempDirectory != null) {
            File fastqDirectory = new File(tempDirectory);
            if (!fastqDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myFastqFileS = DirectoryCrawler.listFileNames
                ("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", fastqDirectory.getAbsolutePath());
            //    (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char'
            // NOTE: If you add addtional file naming conventions here, you must also add them to the "outfile = new File(outFileS.replaceAll()" list below (near the bottom)
            if (myFastqFileS.length == 0 || myFastqFileS == null) {
                printUsage();
                throw new IllegalArgumentException(
                    "Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory: "
                    + tempDirectory);
            } else {
                myLogger.info("FastqToTBTPlugin: setParameters: Using the following fastq files:");
                for (String filename : myFastqFileS) {
                    myLogger.info(filename);
                }
            }
        }
        if (myArgsEngine.getBoolean("-k")) {
            myKeyFile = myArgsEngine.getString("-k");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a key file (option -k).");
        }
        if (myArgsEngine.getBoolean("-e")) {
            myEnzyme = myArgsEngine.getString("-e");
        } else {
            System.out.println("No enzyme specified.  Using enzyme listed in key file.");
//            printUsage();
//            throw new IllegalArgumentException("Please specify the enzyme used to create the GBS library.");
        }
        if (myArgsEngine.getBoolean("-o")) {
            myOutputDir = myArgsEngine.getString("-o");
            File outDirectory = new File(myOutputDir);
            if (!outDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output name you supplied (option -o) is not a directory: " + myOutputDir);
            }
            outDirectory = null;
        }else{
            printUsage();
            throw new IllegalArgumentException("Please specify an output directory (option -o).");
        }
        if (myArgsEngine.getBoolean("-c")) {
            myMinCount = Integer.parseInt(myArgsEngine.getString("-c"));
        } else {
            myMinCount = 1;
        }
        if (myArgsEngine.getBoolean("-y")) { useTBTByte = true; }
        if (myArgsEngine.getBoolean("-y")) {
            if (myArgsEngine.getBoolean("-sh")) {
                printUsage();
                throw new IllegalArgumentException("Options -y and -sh are mutually exclusive.");
            }
            useTBTByte = true;
        } else if (myArgsEngine.getBoolean("-sh")) {
            if (myArgsEngine.getBoolean("-y")) {
                printUsage();
                throw new IllegalArgumentException("Options -y and -sh are mutually exclusive.");
            }
            useTBTShort = true;
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a tagCounts file (-t) *OR* a TagsOnPhysicalMap file (-m)");
        }

        // Create Tags object from tag count file with option -t, or from TOPM file with option -m
        if (myArgsEngine.getBoolean("-t")) {
            if (myArgsEngine.getBoolean("-m")) {
                printUsage();
                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
            }
            myMasterTags = new TagCounts(myArgsEngine.getString("-t"), FilePacking.Bit);
        } else if (myArgsEngine.getBoolean("-m")) {
            if (myArgsEngine.getBoolean("-t")) {
                printUsage();
                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
            }
            myMasterTags = new TagsOnPhysicalMap(myArgsEngine.getString("-m"), true);
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a tagCounts file (-t) *OR* a TagsOnPhysicalMap file (-m)");
        }
    }


    /**
     * Uses an existing Tags object to create one TagsByTaxa file for each fastq file in the input directory.
     *
     * Output TBT files written to the outputDir, using fastq file names with extension changed to .tbt.bin (or .tbt.txt)
     *
     * @param fastqFileS      Array of fastq file names (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param keyFileS       A key file (list of taxa by barcode, lane & flow cell, including plate maps)
     * @param enzyme         The enzyme used to make the library (currently ApeKI or PstI)
     * @param theMasterTags  A Tags object: list of tags to be included in the final TBT
     * @param outputDir      String containing the path of the output directory to contain tags-by-taxa files
     * @param minCount       The minimum number of times a tag must show up in a fastq file before it is included in the corresponding TBT file
     */
    public static void matchTagsToTaxa(String[] fastqFileS, String keyFileS, String enzyme, Tags theMasterTags, String outputDir, int minCount, boolean useTBTByte, boolean useTBTShort) {
        for(int laneNum=0; laneNum<fastqFileS.length; laneNum++) {
            
            //Determine name of output file based on input parameters
            File outfile;
            FilePacking outFormat = useTBTByte ? FilePacking.Byte : (useTBTShort ? FilePacking.Short : FilePacking.Bit);
            String outFileS = outputDir + fastqFileS[laneNum].substring(fastqFileS[laneNum].lastIndexOf(File.separator));
            String replaceS = (outFormat == FilePacking.Text) ? ".tbt.txt" : ((outFormat == FilePacking.Byte) ? ".tbt.byte" : ((outFormat == FilePacking.Short) ? ".tbt.shrt" : ".tbt.bin"));
            outfile = new File(outFileS.replaceAll(
                "(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", 
                replaceS));

            //Skip input file if a corresponding output file has already been written.
            if(outfile.isFile()){
                System.out.println(
                        "An output file "+outfile.getName()+"\n"+ 
                        " already exists in the output directory for file "+fastqFileS[laneNum]+".  Skipping.");
                continue;
            }

            System.out.println("\nWorking on fastq file: " + fastqFileS[laneNum]);
            TagsByTaxa theTBT=null;
            System.gc();
            int goodBarcodedReads=0, allReads=0, goodMatched=0;
            File fastqFile=new File(fastqFileS[laneNum]);
            String[] np=fastqFile.getName().split("_");

            //Create a new object to hold barcoded tags.  The constructor can optionally process a group of fastq
            //files.  A minimum quality score for inclusion of a read can also be provided.
            ParseBarcodeRead thePBR;
            if(np.length==3) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);}
            else if(np.length==5) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);}
            else if(np.length==4) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[0], np[2]);}
            else if(np.length==6) {thePBR=new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);}
            else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_fastq.txt OR code_flowcell_s_lane_fastq.txt");
                System.out.println("   Filename: "+fastqFileS[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
            if(thePBR.getBarCodeCount() == 0){
                System.out.println("No barcodes found.  Skipping this flowcell."); continue;
            }

            //Fill an array with taxon names.
            String[] taxaNames=new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();
            }
            System.out.println("Taxa Names: " + Arrays.toString(taxaNames));

            if (useTBTByte) {
                theTBT=new TagsByTaxaByte(taxaNames,theMasterTags);
            } else if (useTBTShort) {
                theTBT=new TagsByTaxaShort(taxaNames,theMasterTags);
            } else {
                theTBT=new TagsByTaxaBit(taxaNames,theMasterTags);
            }

            // Read the fastq file and assign reads to tags and taxa
            String temp="";
            int currLine = 0;
            goodBarcodedReads=0; allReads=0; goodMatched=0;
            int[]taxaCount= new int[taxaNames.length]; //added by KLS to hold good barcoded read counts at the index of their associated taxa in taxaNames[]
            int[]taxaCountMatch= new int[taxaNames.length];//added by KLS to hold good barcoded reads matched to a tag kept in FastqToTagCount MergeTagCounts
            try{
                BufferedReader br;
                //Read in fastq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if(fastqFileS[laneNum].endsWith(".gz")){
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFileS[laneNum]))));
                }else{
                    br=new BufferedReader(new FileReader(fastqFileS[laneNum]),65536);
                }
                String sl="", qualS="";
                while (((temp = br.readLine()) != null)&&(goodBarcodedReads< maxGoodReads)) {
                    currLine++;
                    //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
                    if((currLine+2)%4==0) {
                        sl = temp;
                    } else if (currLine%4==0) {
                        qualS = temp;
                        allReads++;
                        //After quality score is read, decode barcode using the current sequence & quality  score
                        ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sl, qualS, true, 0);
                        if(rr!=null) {
                            goodBarcodedReads++;
                            int t=theTBT.getIndexOfTaxaName(rr.getTaxonName());
                            int h=theTBT.getTagIndex(rr.getRead());
                            taxaCount[t]++;//add one to index of taxa; counts goodBarcodedReads
                            if(h>-1) {
                                theTBT.addReadsToTagTaxon(h, t, 1);
                                taxaCountMatch[t]++; //KLS add one to index of taxa; counts good barcoded reads matched to a tag based on count mimimum specified in tagCounts
                                goodMatched++;
                            }
                        }
                        if(allReads%1000000==0) System.out.println("Total Reads:"+allReads+" goodReads:"+goodBarcodedReads+" goodMatched:"+goodMatched);
                        if(allReads%1000000==0) System.out.println("#goodBarcodes: " + Arrays.toString(taxaCount)); //KLS prints out goodbarcoded read counts in same order as taxa names
                    }
                 }
                br.close();
            }
            catch(Exception e) {
                System.out.println("Catch testBasicPipeline c="+goodBarcodedReads+" e="+e);
                System.out.println(temp);
                e.printStackTrace();
            }
            System.out.println("Timing process (writing TagsByTaxa file)..."); long timePoint1 = System.currentTimeMillis();
            theTBT.writeDistFile(outfile,outFormat, minCount);
            System.out.println("...process (writing TagsByTaxa file) took "+(System.currentTimeMillis()-timePoint1)+" milliseconds.");
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Total number of good, barcoded reads="+goodBarcodedReads);
            System.out.println("Taxa: " + Arrays.toString(taxaNames)); //KLS prints out  taxa names
            System.out.println("#goodBarcodes: " + Arrays.toString(taxaCount)); //KLS prints out goodbarcoded read counts in same order as taxa names
            int filesDone = laneNum + 1;
            System.out.println("Finished reading "+filesDone+" of "+fastqFileS.length+" sequence files: "+fastqFileS[laneNum]+"\n");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
