/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.github.codingtogenomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import net.minidev.json.parser.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 *
 * @author david
 */
public class CodingToGenomic {

    /**
     * @param args the command line arguments
     */
  public static String SERVER = "http://rest.ensembl.org";
  public static final String GRCh37Server = "http://grch37.rest.ensembl.org/";
  public static final JSONParser PARSER = new JSONParser(JSONParser.MODE_JSON_SIMPLE);

  public static int requestCount = 0;
  public static long lastRequestTime = System.currentTimeMillis();

  public static void main(String[] args) throws Exception {
      
      //parse commandline
    Options options = new Options();
    CommandLineParser parser = new PosixParser();
    String gene = new String();
    String transcript = new String();
    String species = "human";
    String coordinate = new String();
    StringBuilder errorMsg = new StringBuilder();
    try{
        options = getOptions(args); 
    }catch(org.apache.commons.cli.ParseException ex){
        System.err.println( "Parsing failed.  Reason: " + ex.getMessage() );
        System.exit(1);
    }
    CommandLine line = parser.parse(options, args);
    if (line.hasOption("help")){
        showHelp(options);
    }
    if (line.hasOption("gene")){
        gene = line.getOptionValue("gene");
    }else{
        if (! line.hasOption("transcript")){
            errorMsg.append("Either --gene or --transcript argument is required\n");
        }
    }
    if (line.hasOption("transcript")){
        if (line.hasOption("gene")){
            errorMsg.append("Please specify only one of "
                    + "--gene or --transcript arguments, not both\n");
        }else{
            transcript = line.getOptionValue("transcript");
            if (line.hasOption("species")){
                System.out.println("Ignoring --species option when using --transcript argument");
            }
        }
    }
    if (line.hasOption("coordinate")){
        coordinate = line.getOptionValue("coordinate");
    }else{
        errorMsg.append("--coordinate argument is required\n");
    }
    if (line.hasOption("species")){
        species = line.getOptionValue("species");
    }
    if (line.hasOption("b37")){
        if (species.equalsIgnoreCase("human") || species.equalsIgnoreCase("homo sapiens")){
            SERVER = GRCh37Server;
        }else{
            System.out.println("--b37 argument will be ignored - it can only be "
                    + "used when human is the species of interest. Current species"
                    + " is " + species + ".\n");
        }
    }
    
    if (errorMsg.length() > 0){
        showHelp(options, errorMsg.toString(), 2);
    }
    int c = 0;
    try{
        c = Integer.parseInt(coordinate);
    }catch (NumberFormatException ex){
        showHelp(options, "--coordinate argument '" + coordinate + "' is not an integer", 2);
    }

    //got arguments
    if (! gene.isEmpty()){
        codingToGenomic(species, gene, c);
    }else{
        codingToGenomicTranscript(species, transcript, c);
    }
    
    
    /*ENSEMBL EXAMPLE
    
    JSONArray variants = getVariants(species, gene);
    for(Object variantObject: variants) {
      JSONObject variant = (JSONObject)variantObject;
      String srName = (String)variant.get("seq_region_name");
      Number start = (Number)variant.get("start");
      Number end = (Number)variant.get("end");
      Number strand = (Number)variant.get("strand");
      String id = (String)variant.get("id");
      String consequence = (String)variant.get("id");
      String output = String.format("%s:%d-%d:%d ==> %s (%s)", srName, start, end, strand, id, consequence);
      System.out.println(output);
    }
            */
  }
  
  public static List<String> getGeneAndSymbolFromTranscript(String id)throws ParseException, MalformedURLException, IOException, InterruptedException {
      String endpoint = "/overlap/id/ENST00000398208?feature=gene";
      JSONArray genes = (JSONArray) getJSON(endpoint);
      JSONObject gene = (JSONObject)genes.get(0);
      String ensid = (String) gene.get("id");
      String name = (String) gene.get("external_name");
      return Arrays.asList(ensid, name);
  }
  
  public static void codingToGenomicTranscript(String species, String id, int c) throws ParseException, MalformedURLException, IOException, InterruptedException {
      String endpoint = "/lookup/id/"+id+"?expand=1";
      JSONObject tr = (JSONObject) getJSON(endpoint);
      if (tr.isEmpty()){
        System.out.println("No protein coding transcripts found for " + id);
        return;
    }
    List<String> geneAndSymbol = getGeneAndSymbolFromTranscript(id);  
      
    String seq = getCds(id);
    if (seq != null){
        if (seq.length() >= c ){
            String gCoord = cdsToGenomicCoordinate(id, c);
            if (gCoord != null){
                System.out.println(geneAndSymbol.get(1) +" " + geneAndSymbol.get(0) 
                        + " " + id + " c." + c + " => " + gCoord);
            }
        }else{
            System.out.println("CDS coordinate " + c + " is greater than "
                    + "length of CDS (" + seq.length() + ") for " + tr);
        }
    }else{
        System.out.println("ERROR: No CDS sequence found for " + tr);
    }
  }
  
  public static void codingToGenomic(String species, String symbol, int c) throws ParseException, MalformedURLException, IOException, InterruptedException {
    String id = getGeneID(species, symbol);
    ArrayList<String> tr = getTranscriptIds(id);
    if (tr.isEmpty()){
        System.out.println("No protein coding transcripts found for " + symbol
                + " (" + id + ")");
        return;
    }
//    System.out.println(id + " => " + String.join(",", tr));
    for (String t : tr){
        String seq = getCds(t);
        if (seq != null){
            if (seq.length() >= c ){
                String gCoord = cdsToGenomicCoordinate(t, c);
                if (gCoord != null){
                    System.out.println(symbol + " " + id + " " + t + 
                            " c." + c + " => " + gCoord);
                }
            }else{
                System.out.println("CDS coordinate " + c + " is greater than "
                        + "length of CDS (" + seq.length() + ") for " + t);
            }
        }else{
            System.out.println("ERROR: No CDS sequence found for " + t);
        }
    }
  }
  
  public static String getCds(String id) throws ParseException, MalformedURLException, IOException, InterruptedException {
      String endpoint = "/sequence/id/" + id + "?type=cds";
      JSONObject sequence = (JSONObject) getJSON(endpoint);
      if (sequence.containsKey("seq")){
          return (String) sequence.get("seq");
      }
      return null;
  }
  
  public static String cdsToGenomicCoordinate(String id, int coord) throws ParseException, MalformedURLException, IOException, InterruptedException {
      String endpoint = "/map/cds/" + id +"/"+ coord + ".." + coord;
      JSONObject info = (JSONObject) getJSON(endpoint);
      StringBuilder mapStrings = new StringBuilder();
      if(info.isEmpty()) {
          throw new RuntimeException("Got nothing for endpoint "+endpoint);
      }
      if (info.containsKey("mappings")){
          JSONArray mappings = (JSONArray) info.get("mappings");
          for (Object map: mappings){
              JSONObject m = (JSONObject) map;
              if (m.containsKey("start") && m.containsKey("seq_region_name")){
                  String assembly = (String) m.get("assembly_name");
                  String chr = (String) m.get("seq_region_name");
                  Long start = (Long) m.get("start");
                  mapStrings.append(chr + ":" + start + " (" + assembly + ")\n");
              }
          }
          return mapStrings.toString();
      }
      return null;
  }
  
  public static ArrayList<String> getTranscriptIds(String id) throws ParseException, MalformedURLException, IOException, InterruptedException {
    ArrayList<String> transcriptIds = new ArrayList<>();
    String endpoint = "/lookup/id/"+id+"?expand=1";
    JSONObject info = (JSONObject) getJSON(endpoint);
    if(info.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    if (info.containsKey("Transcript")){
        JSONArray trs = (JSONArray) info.get("Transcript");
        for (Object t: trs){
           JSONObject j = (JSONObject) t;
           String biotype = (String)j.get("biotype");
           if (biotype.equals("protein_coding")){
               transcriptIds.add((String) j.get("id"));
           }
        }
    }
    return transcriptIds;
  }

  
  public static JSONArray getVariants(String species, String symbol) throws ParseException, MalformedURLException, IOException, InterruptedException {
    String id = getGeneID(species, symbol);
    return (JSONArray) getJSON("/overlap/id/"+id+"?feature=variation");
  }

  public static String getGeneID(String species, String symbol) throws ParseException, MalformedURLException, IOException, InterruptedException {
    String endpoint = "/xrefs/symbol/"+species+"/"+symbol+"?object_type=gene";
    JSONArray genes = (JSONArray) getJSON(endpoint);
    if(genes.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    JSONObject gene = (JSONObject)genes.get(0);
    return (String)gene.get("id");
  }

  public static Object getJSON(String endpoint) throws ParseException, MalformedURLException, IOException, InterruptedException {
    String jsonString = getContent(endpoint);
    return PARSER.parse(jsonString);
  }

  public static String getContent(String endpoint) throws MalformedURLException, IOException, InterruptedException {

    if(requestCount == 15) { // check every 15
      long currentTime = System.currentTimeMillis();
      long diff = currentTime - lastRequestTime;
      //if less than a second then sleep for the remainder of the second
      if(diff < 1000) {
        Thread.sleep(1000 - diff);
      }
      //reset
      lastRequestTime = System.currentTimeMillis();
      requestCount = 0;
    }

    URL url = new URL(SERVER+endpoint);
    URLConnection connection = url.openConnection();
    HttpURLConnection httpConnection = (HttpURLConnection)connection;
    httpConnection.setRequestProperty("Content-Type", "application/json");

    InputStream response = httpConnection.getInputStream();
    int responseCode = httpConnection.getResponseCode();

    if(responseCode != 200) {
      if(responseCode == 429 && httpConnection.getHeaderField("Retry-After") != null) {
        double sleepFloatingPoint = Double.valueOf(httpConnection.getHeaderField("Retry-After"));
        double sleepMillis = 1000 * sleepFloatingPoint;
        Thread.sleep((long)sleepMillis);
        return getContent(endpoint);
      }
      throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
    }

    String output;
    Reader reader = null;
    try {
      reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      StringBuilder builder = new StringBuilder();
      char[] buffer = new char[8192];
      int read;
      while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
        builder.append(buffer, 0, read);
      }
      output = builder.toString();
    } 
    finally {
      if (reader != null) {
        try {
          reader.close(); 
        } 
        catch (IOException logOrIgnore) {
          logOrIgnore.printStackTrace();
        }
      }
    }

    return output;
  }
  
  private static Options getOptions (String args[]) throws org.apache.commons.cli.ParseException{
      Option help = new Option( "h", "help", false, "print this message" );
      Option b37 = new Option( "b", "b37", false,  "use b37/hg19 coordinates (human only)" );
      Option gene   = OptionBuilder.withLongOpt("gene")
                                .withArgName( "gene" )
                                .hasArg()
                                .withDescription(  "Gene symbol to search" )
                                .create( "g" );
      Option transcript   = OptionBuilder.withLongOpt("transcript")
                                .withArgName( "gene" )
                                .hasArg()
                                .withDescription(  "Ensembl transcript ID to search" )
                                .create( "t" );
      Option species   = OptionBuilder.withLongOpt("species")
                                .withArgName( "species" )
                                .hasArg()
                                .withDescription(  "Name of species to search"
                                        + " (default = human)" )
                                .create( "s" );
      Option coordinate   = OptionBuilder.withLongOpt("coordinate")
                                .withArgName( "coordinate" )
                                .hasArg()
                                .withDescription(  "CDS coordinate" )
                                .create( "c" );
      Options options = new Options();
      options.addOption(gene);
      options.addOption(species);
      options.addOption(coordinate);
      options.addOption(b37);
      options.addOption(transcript);
      options.addOption(help);
      return options;
  }
  
  private static void showHelp(Options o){
      showHelp(o, null, 0);
  }
  
  private static void showHelp(Options  o, String msg, Integer exitVal){
      HelpFormatter formatter = new HelpFormatter();
      formatter.setOptionComparator(new OptionComparator());
      formatter.setLeftPadding(8);
      formatter.setWidth(80);
      if (msg != null){
          System.out.println("\nERROR: " + msg);
      }
      formatter.printHelp("java -jar CodingToGenomic.jar "
              + "[-g <gene>] [-c <coordinate>] [options]\n", o, false);
      System.exit(exitVal);
  }
  
  static class OptionComparator<T extends Option> implements Comparator<T> {

        private static final List<String> OPT_ORDER = Arrays.asList(
                "g", "t", "c", "s", "b", "h");

        public int compare(T o1, T o2) {
            return OPT_ORDER.indexOf(o1.getOpt()) - OPT_ORDER.indexOf(o2.getOpt());
        }
    }
}
