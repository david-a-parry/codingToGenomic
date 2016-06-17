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
import java.util.Collections;
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
  private static String SERVER = "http://rest.ensembl.org";
  private static final String GRCh37Server = "http://grch37.rest.ensembl.org/";
  private static final JSONParser PARSER = new JSONParser(JSONParser.MODE_JSON_SIMPLE);

  private static int requestCount = 0;
  private static long lastRequestTime = System.currentTimeMillis();

  public static void main(String[] args) throws Exception {
      
      //parse commandline
    Options options = new Options();
    CommandLineParser parser = new PosixParser();
    String gene = new String();
    String transcript = new String();
    String species = "human";
    boolean mapCdna = false;
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
        species = line.getOptionValue("species").replaceAll("\\s+", "_");
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
    if (line.hasOption("noncoding")){
        mapCdna = true;
    }
    
    if (errorMsg.length() > 0){
        showHelp(options, errorMsg.toString(), 2);
    }
    int c = 0;
    boolean threePrimeUtr = false;
    String prefix = "c.";
    if (mapCdna){
        prefix = "n.";
        try{
            c = Integer.parseInt(coordinate);
        }catch (NumberFormatException ex){
            showHelp(options, "--coordinate argument '" + coordinate + "' could not "
                    + "be parsed as an integer", 2);
        }
    }else if (coordinate.startsWith("*")){
        threePrimeUtr = true;
        prefix = "c.*";
        String coord = coordinate.replaceFirst("\\*", "");
        try{
            c = Integer.parseInt(coord);
        }catch (NumberFormatException ex){
            showHelp(options, "--coordinate argument '" + coordinate + "' could not "
                    + "be parsed as an integer or UTR coordinate", 2);
        }
    }else{
        try{
            c = Integer.parseInt(coordinate);
        }catch (NumberFormatException ex){
            showHelp(options, "--coordinate argument '" + coordinate + "' could not "
                    + "be parsed as an integer", 2);
        }
    }
    //got arguments
    String result;
    String header = "Input\tSymbol\tEnsemblGene\tEnsemblTranscript\tGenomicCoordinate";
    if (! gene.isEmpty()){
        IdParser idParser = new IdParser(gene);
        System.out.println("Interpretting " + gene + " as of type " + idParser.getIdentifierType());
        if (idParser.isEnsemblId()){
            if (line.hasOption("species")){
                System.out.println("Ignoring --species option when searching Ensembl ID.");
            }
            if (idParser.isTranscript()){
                result = codingToGenomicTranscript(gene, c, threePrimeUtr, mapCdna);
            }else if (idParser.isEnsp()){
                result = codingToGenomicEnsp(gene, c,threePrimeUtr, mapCdna);
            }else{
                result = codingToGenomicId(gene, c,threePrimeUtr, mapCdna);
            }
        }else{
            if (idParser.isTranscript()){
                //append user input to beginning
                result = codingToGenomicXrefTranscript(species, gene, c, 
                        threePrimeUtr, mapCdna);
            }else{
                result = codingToGenomicXref(species, gene, c, 
                        threePrimeUtr, mapCdna);
            }
        }
        if (idParser.isTranscript() || idParser.isEnsp()){
            
            result = gene + ":" + prefix + c + "\t" + result;
        }else{
            result = convertGeneResult(result, gene, c, prefix);
        }
        
    }else{
        System.out.println("Searching for " + transcript + " as Ensembl transcript ID") ; 
        result = codingToGenomicTranscript(transcript, c, threePrimeUtr, mapCdna);
        //append user input to beginning
        result = transcript + ":" + prefix + c + "\t" + result;
    }
    
    System.out.println(header);
    System.out.println(result);
    
  }
  
  private static String convertGeneResult(final String r, final String input, 
          final int c, final String prefix){
    String split[] = r.split("\n");
    StringBuilder converted = new StringBuilder();
    for (String s: split){
        converted.append(input).append(":").append(prefix).append(c).append("\t").append(s).append("\n");
    }
    return converted.toString();
  }
  
  private static List<String> getGeneAndSymbolFromTranscript(final String id)throws ParseException, MalformedURLException, IOException, InterruptedException {
        String endpoint = "/lookup/id/" + id ;
        JSONObject gene = (JSONObject) getJSON(endpoint);
        String symbol = new String();
        if (gene.containsKey("display_name")){
            String display_name = (String) gene.get("display_name");
            symbol = display_name.split("-")[0];
        }
        endpoint = "/overlap/id/" + id + "?feature=gene";
        JSONArray genes = (JSONArray) getJSON(endpoint);
        if (genes.isEmpty()){
            throw new RuntimeException(String.format("Could not get gene details "
                  + "for transcript %s.\nRetrieval from URL %s%s returned nothing.", 
                  id, SERVER, endpoint));
        }
        for (Object o: genes){
            JSONObject g = (JSONObject) o;
            String name = (String) g.get("external_name");
            if (name.equalsIgnoreCase(symbol)){
                String ensid = (String) g.get("id");
                return (Arrays.asList(ensid, symbol));
            }
        }
        //failed to get ensgene id 
        return (Arrays.asList("", symbol));
    }
  
  private static String codingToGenomicTranscript(final String id, final int c, 
          final boolean utr, final boolean mapCdna) throws ParseException,
          MalformedURLException, IOException, InterruptedException {
    String endpoint = "/lookup/id/"+id+"?expand=1";
    JSONObject tr = (JSONObject) getJSON(endpoint);
    if (tr.isEmpty()){
        return "-\t-\t-\t" + id + "\tNo transcripts found for " + id;
    }
    List<String> geneAndSymbol = getGeneAndSymbolFromTranscript(id);  
    String stub = geneAndSymbol.get(1) +"\t" + geneAndSymbol.get(0) 
                        + "\t" + id + "\t" ;
    String fetch = "cds";
    if (mapCdna){
        fetch = "cdna";
    }
    String seq = getTranscriptSequence(id, fetch);
    if (seq != null){
        if (seq.length() >= c ){
            String gCoord;
            if (mapCdna){
                gCoord = cdnaToGenomicCoordinate(id, c);
            }else{
                gCoord = cdsToGenomicCoordinate(id, c, utr);
            }
            if (gCoord != null){
                return stub + gCoord + "\n";
            }
        }else{
            return stub + fetch.toUpperCase() + " coordinate (" + c + ") > " + 
                    fetch.toUpperCase() + "length (" + seq.length() + ")\n";
        }
    }
    return stub + "No " + fetch.toUpperCase() + " sequence found\n";
  }
  
  private static String codingToGenomicId(final String id, final int c, 
          final boolean utr, final boolean mapCdna)throws ParseException, MalformedURLException, 
          IOException, InterruptedException {
      final String symbol = getGeneSymbol(id);
      return codingToGenomicGene(id, symbol, c, utr, mapCdna);
  }
  
  private static String codingToGenomicEnsp(final String id, final int c, 
          final boolean utr, final boolean mapCdna)throws ParseException, MalformedURLException, 
          IOException, InterruptedException {
      //get parent Transcript from ENSP ID and process as transcript...
      final String transcript = getTranscriptFromEnsp(id);
      if (transcript == null){
          return "-\t-\t-\t-Could not identify parent transcript\n";
      }else{
          return codingToGenomicTranscript(transcript, c, utr, mapCdna); 
      }
  }
  
  
  private static String codingToGenomicXref(final String species, final String id, 
          final int c, final boolean utr, final boolean mapCdna) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
      final String ensid = getGeneID(species, id);
      final String symbol = getGeneSymbol(ensid);
      return codingToGenomicGene(ensid, symbol, c, utr, mapCdna);
  }
  
  private static String codingToGenomicXrefTranscript(final String species, 
          final String id, final int c, final boolean utr, final boolean mapCdna) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
      final String ensid = getTranscriptXrefId(species, id);
      //String symbol = getGeneSymbol(ensid);
      return codingToGenomicTranscript(ensid, c, utr, mapCdna);
  }
  
  private static String codingToGenomicSymbol(final String species, 
          final String symbol, final int c, final boolean utr, final boolean mapCdna) 
          throws ParseException, MalformedURLException, IOException, 
          InterruptedException {
      final String id = getGeneID(species, symbol);
      return codingToGenomicGene(id, symbol, c, utr, mapCdna);
  }
  
  private static String codingToGenomicGene(final String id, final String symbol, 
          final int c, final boolean utr, final boolean mapCdna) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
    final ArrayList<String> tr = getTranscriptIds(id, mapCdna);
    if (tr.isEmpty()){
        return "-\t-\t-\t" + id + "\tNo coding transcripts found for " + symbol
                + " (" + id + ")\n";
    }
    /*
    StringBuilder transcriptList = new StringBuilder(tr.get(0));
    for (int i = 1; i < tr.size(); i++){
        transcriptList.append(",").append(tr.get(i));
    }
    //System.out.println(symbol + " => " + id + " => " + transcriptList.toString());
    */
    String fetch = "cds";
    if (mapCdna){
        fetch = "cdna";
    }
    final String stub = symbol + "\t" + id + "\t";
    StringBuilder results = new StringBuilder();
    for (String t : tr){
        final String seq = getTranscriptSequence(t, fetch);
        if (seq != null){
            if (seq.length() >= c ){
                String gCoord;
                if (mapCdna){
                    gCoord = cdnaToGenomicCoordinate(t, c);
                }else{
                    gCoord = cdsToGenomicCoordinate(t, c, utr);
                }
                if (gCoord != null){
                    results.append(stub).append(t).append("\t").append(gCoord)
                            .append("\n");
                }else{
                    results.append(stub).append(t).append("\tCould not map ")
                            .append(fetch.toUpperCase()).append(" coordinate\n");
                }
            }else{
                results.append(stub).append(t).append("\t")
                        .append(fetch.toUpperCase()).append("coordinate (")
                        .append(c).append(") > ").append(fetch.toUpperCase())
                        .append(" length (").append(seq.length()).append(")\n");
            }
        }else{
            results.append(stub).append("-\tNo CDS sequence found\n");
        }
        
    }
    return results.toString();
  }
  
  private static String getTranscriptFromEnsp(final String id) throws 
          ParseException, MalformedURLException, IOException, InterruptedException {
    final String endpoint = "/lookup/id/"+id;
    final JSONObject info = (JSONObject) getJSON(endpoint);
    if(info.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    if (info.containsKey("Parent")){
        return (String) info.get("Parent");
    }else{
        return null;
    }
  }
  
  
  private static String getTranscriptSequence(final String id, final String type)
          throws ParseException, MalformedURLException, IOException, InterruptedException {
      final String endpoint = "/sequence/id/" + id + "?type=" + type;
      final JSONObject sequence = (JSONObject) getJSON(endpoint);
      if (sequence.containsKey("seq")){
          return (String) sequence.get("seq");
      }
      return null;
  }
  
  private static String cdsToGenomicCoordinate(final String id, final int coord, 
          final boolean utr) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
      
      if (coord < 0 || utr){
          return utrToGenomicCoordinate(id, coord);
      }
      
      final String endpoint = "/map/cds/" + id +"/"+ coord + ".." + coord;
      final JSONObject info = (JSONObject) getJSON(endpoint);
      StringBuilder mapStrings = new StringBuilder();
      if(info.isEmpty()) {
          throw new RuntimeException("Got nothing for endpoint "+endpoint);
      }
      if (info.containsKey("mappings")){
          final JSONArray mappings = (JSONArray) info.get("mappings");
          int n = 0;
          for (Object map: mappings){
              n++;
              if (n > 1){
                  mapStrings.append("|");
              }
              final JSONObject m = (JSONObject) map;
              if (m.containsKey("start") && m.containsKey("seq_region_name")){
                  final String assembly = (String) m.get("assembly_name");
                  final String chr = (String) m.get("seq_region_name");
                  final Long start = (Long) m.get("start");
                  mapStrings.append(chr).append(":").append(start).append(" (").append(assembly).append(")");
              }
          }
          return mapStrings.toString();
      }
      return null;
  }
  
  private static String utrToGenomicCoordinate(final String id, final int coord) 
          throws ParseException, MalformedURLException, IOException, InterruptedException {
      if (coord > 0){
          return threePrimeUtrToGenomicCoordinate(id, coord);
      }else{
          return fivePrimeUtrToGenomicCoordinate(id, coord);
      }
  }
  
  private static String threePrimeUtrToGenomicCoordinate(final String id, final 
          int coord) throws ParseException, MalformedURLException, IOException, 
          InterruptedException {
      final String trEndpoint = "/lookup/id/"+id+"?expand=1";
      final JSONObject trInfo = (JSONObject) getJSON(trEndpoint);
      if(trInfo.isEmpty()) {
        throw new RuntimeException("Got nothing for endpoint "+trEndpoint);
      }
      final TranscriptDetails trDetails = getTranscriptDetailsFromJson(trInfo);
      final int cdsLength = trDetails.getCodingLength();
      final String cdsEndpoint = "/map/cds/" + id +"/"+ cdsLength + ".." + cdsLength;
      final JSONObject cdsInfo = (JSONObject) getJSON(cdsEndpoint);

      if(cdsInfo.isEmpty()) {
        throw new RuntimeException("Got nothing for endpoint "+cdsEndpoint);
      }
      
      if (cdsInfo.containsKey("mappings")){
          StringBuilder mapStrings = new StringBuilder();
          final JSONArray mappings = (JSONArray) cdsInfo.get("mappings");
          int n = 0;
          for (Object map: mappings){
              n++;
              if (n > 1){
                  mapStrings.append("|");
              }
              JSONObject m = (JSONObject) map;
              if (m.containsKey("start") && m.containsKey("seq_region_name")){
                  //assembly = (String) m.get("assembly_name");
                  String chr = (String) m.get("seq_region_name");
                  Long genomicCoord = (Long) m.get("start");
                  String trPos = trDetails.getCdnaPosition(chr, genomicCoord.intValue());
                  try{
                      int pos = Integer.parseInt(trPos);
                      if (coord > pos ){
                          mapStrings.append("UTR position > length of 3' UTR");
                          
                      }else{
                          String cdnaPos = cdnaToGenomicCoordinate(id, pos + coord);
                          if(cdnaPos != null){
                              mapStrings.append(cdnaPos);
                          }
                      }
                  }catch(NumberFormatException | ParseException | IOException | InterruptedException ex){
                      mapStrings.append( "Could not parse cDNA position "
                              + "(" + trPos + ") for " + id);
                  }
              }
          }
          return mapStrings.toString();
      }
      return null;
  }
  
  private static String fivePrimeUtrToGenomicCoordinate(final String id, final int coord) throws ParseException, MalformedURLException, IOException, InterruptedException {
  
      if (coord > 0){
          throw new RuntimeException("fivePrimeUtrToGenomicCoordinate method requires " +
                  "coordinate to be less than zero");
      }

      final String trEndpoint = "/lookup/id/"+id+"?expand=1";
      final JSONObject trInfo = (JSONObject) getJSON(trEndpoint);
      if(trInfo.isEmpty()) {
        throw new RuntimeException("Got nothing for endpoint "+trEndpoint);
      }
      final TranscriptDetails trDetails = getTranscriptDetailsFromJson(trInfo);

      final String cdsEndpoint = "/map/cds/" + id +"/"+ 1 + ".." + 1;
      final JSONObject cdsInfo = (JSONObject) getJSON(cdsEndpoint);

      if(cdsInfo.isEmpty()) {
        throw new RuntimeException("Got nothing for endpoint "+cdsEndpoint);
      }
      
      if (cdsInfo.containsKey("mappings")){
          StringBuilder mapStrings = new StringBuilder();
          final JSONArray mappings = (JSONArray) cdsInfo.get("mappings");
          int n = 0;
          for (Object map: mappings){
              n++;
              if (n > 1){
                  mapStrings.append("|");
              }
              JSONObject m = (JSONObject) map;
              if (m.containsKey("start") && m.containsKey("seq_region_name")){
                  //assembly = (String) m.get("assembly_name");
                  String chr = (String) m.get("seq_region_name");
                  Long genomicCoord = (Long) m.get("start");
                  String trPos = trDetails.getCdnaPosition(chr, genomicCoord.intValue());
                  try{
                      int pos = Integer.parseInt(trPos);
                      if (Math.abs(coord) > pos ){
                          mapStrings.append("UTR position > length of 5' UTR");
                      }else{
                          String cdnaPos = cdnaToGenomicCoordinate(id, pos + coord);
                          if(cdnaPos != null){
                              mapStrings.append(cdnaPos);
                          }
                      }
                  }catch(NumberFormatException | ParseException | IOException | InterruptedException ex){
                      mapStrings.append( "Could not parse cDNA position "
                              + "(" + trPos + ") for " + id);
                  }
              }
          }
          return mapStrings.toString();
      }
      return null;
  }
  
  
    private static String cdnaToGenomicCoordinate(final String id, 
            final int coord) throws ParseException, MalformedURLException, 
            IOException, InterruptedException {
      
      if (coord < 0){
          return utrToGenomicCoordinate(id, coord);
      }
      
      final String endpoint = "/map/cdna/" + id +"/"+ coord + ".." + coord;
      final JSONObject info = (JSONObject) getJSON(endpoint);
      StringBuilder mapStrings = new StringBuilder();
      if(info.isEmpty()) {
          throw new RuntimeException("Got nothing for endpoint "+endpoint);
      }
      if (info.containsKey("mappings")){
          final JSONArray mappings = (JSONArray) info.get("mappings");
          int n = 0;
          for (Object map: mappings){
              n++;
              if (n > 1){
                  mapStrings.append("|");
              }
              JSONObject m = (JSONObject) map;
              if (m.containsKey("start") && m.containsKey("seq_region_name")){
                  String assembly = (String) m.get("assembly_name");
                  String chr = (String) m.get("seq_region_name");
                  Long start = (Long) m.get("start");
                  if (start == coord){
                      mapStrings.append("Not found");
                  }else{
                      mapStrings.append(chr + ":" + start + " (" + assembly + ")");                      
                  }
              }
          }
          return mapStrings.toString();
      }
      return null;
    }
  
    private static TranscriptDetails getTranscriptDetailsFromJson(final JSONObject j)
            throws ParseException, MalformedURLException, IOException, InterruptedException {
        final TranscriptDetails trans = new TranscriptDetails();        
        trans.setTranscriptId((String) j.get("id"));
        final String biotype = (String)j.get("biotype");
        trans.setBiotype(biotype);
        if (j.containsKey("is_canonical")){
            String isCanon = j.get("is_canonical").toString();
            if (Integer.parseInt(isCanon) > 0){
                trans.setIsCanonical(true);
            }else{
                trans.setIsCanonical(false);
            }
        }
        if (j.containsKey("strand")){
           if ((Long) j.get("strand") > 0){
               trans.setStrand(1);
           }else{
               trans.setStrand(-1);
           }
       }
       //get exons
       if (j.containsKey("Exon")){
           JSONArray exons = (JSONArray) j.get("Exon");
           for (Object e: exons){
               JSONObject jxon = (JSONObject) e;
               TranscriptDetails.Exon exon = trans.new Exon();
               Long start = (Long) jxon.get("start");
               Long end = (Long) jxon.get("end");
               exon.setStart(start.intValue());
               exon.setEnd(end.intValue());
               trans.getExons().add(exon);
           }
           //sort and number exons
           Collections.sort(trans.getExons());
           for (int i = 0; i < trans.getExons().size(); i++){
               if (trans.getStrand() < 0){
                   trans.getExons().get(i).setOrder(
                           trans.getExons().size() - i);
               }else{
                   trans.getExons().get(i).setOrder(i+1);
               }
           }
       }
       
       //get chromosome
       if (j.containsKey("seq_region_name")){
            trans.setChromosome((String) j.get("seq_region_name"));
        }

      //get transcription start and end
       if (j.containsKey("start")){
           Long start = (Long)j.get("start");
           trans.setTxStart(start.intValue());
       }
       if (j.containsKey("end")){
           Long end = (Long)j.get("end");
           trans.setTxEnd(end.intValue());
       }
       
       
       //get translation start and end if coding                   
       if (j.containsKey("Translation")){
           JSONObject p = (JSONObject) j.get("Translation");
           trans.setProteinId((String) p.get("id"));
           Long start = (Long) p.get("start");
           Long end = (Long) p.get("end");
           trans.setCdsStart(start.intValue());
           trans.setCdsEnd(end.intValue());
           Long length = (Long) p.get("length");
           trans.setProteinLength(length.intValue());
       }else{
           //if using a transcript id for some reason rest won't return translation
       }
       return trans;
    }  
  
  //requires ensembl gene ID for id
  private static ArrayList<String> getTranscriptIds(final String id, 
          final boolean getNonCoding) 
          throws ParseException, MalformedURLException, IOException, InterruptedException {
    final ArrayList<String> transcriptIds = new ArrayList<>();
    final String endpoint = "/lookup/id/"+id+"?expand=1";
    final JSONObject info = (JSONObject) getJSON(endpoint);
    if(info.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    if (info.containsKey("Transcript")){
        final JSONArray trs = (JSONArray) info.get("Transcript");
        for (Object t: trs){
           final JSONObject j = (JSONObject) t;
           final String biotype = (String)j.get("biotype");
           if (biotype.equals("protein_coding") || getNonCoding){
               transcriptIds.add((String) j.get("id"));
           }
        }
    }
    return transcriptIds;
  }
  
  private static String getTranscriptXrefId(final String species, final String xref) 
          throws ParseException, MalformedURLException, IOException, InterruptedException {
    final String endpoint = "/xrefs/symbol/"+species+"/"+xref+"?object_type=transcript";
    final JSONArray trs = (JSONArray) getJSON(endpoint);
    if(trs.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    final JSONObject tr = (JSONObject)trs.get(0);
    return (String)tr.get("id");
  }
  

  
  private static JSONArray getVariants(final String species, final String symbol) 
          throws ParseException, MalformedURLException, IOException, InterruptedException {
    final String id = getGeneID(species, symbol);
    return (JSONArray) getJSON("/overlap/id/"+id+"?feature=variation");
  }

  private static String getGeneID(final String species, final String symbol) 
          throws ParseException, MalformedURLException, IOException, InterruptedException {
    final String endpoint = "/xrefs/symbol/"+species+"/"+symbol+"?object_type=gene";
    final JSONArray genes = (JSONArray) getJSON(endpoint);
    if (genes.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    final JSONObject gene = (JSONObject)genes.get(0);
    return (String)gene.get("id");
  }
  
  private static String getGeneSymbol(final String id) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
    final String endpoint = "/lookup/id/"+id;
    final JSONObject gene = (JSONObject) getJSON(endpoint);
    if (gene.isEmpty()) {
      throw new RuntimeException("Got nothing for endpoint "+endpoint);
    }
    return (String)gene.get("display_name");
  }

  private static Object getJSON(final String endpoint) throws ParseException, 
          MalformedURLException, IOException, InterruptedException {
    final String jsonString = getContent(endpoint);
    return PARSER.parse(jsonString);
  }

  private static String getContent(final String endpoint) throws 
          MalformedURLException, IOException, InterruptedException {

    if(requestCount == 15) { // check every 15
      final long currentTime = System.currentTimeMillis();
      final long diff = currentTime - lastRequestTime;
      //if less than a second then sleep for the remainder of the second
      if(diff < 1000) {
        Thread.sleep(1000 - diff);
      }
      //reset
      lastRequestTime = System.currentTimeMillis();
      requestCount = 0;
    }

    final URL url = new URL(SERVER+endpoint);
    URLConnection connection = url.openConnection();
    HttpURLConnection httpConnection = (HttpURLConnection)connection;
    httpConnection.setRequestProperty("Content-Type", "application/json");

    final InputStream response = httpConnection.getInputStream();
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
      Option noncoding = new Option( "n", "noncoding", false,  "use transcript "
              + "cDNA coordinates, rather than CDS coordinates" );
      Option gene   = OptionBuilder.withLongOpt("gene")
                                .withArgName( "gene" )
                                .hasArg()
                                .withDescription(  "Gene symbol to search" )
                                .create( "g" );
      Option transcript   = OptionBuilder.withLongOpt("transcript")
                                .withArgName( "Ensembl transcript" )
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
                                .withDescription(  "CDS coordinate (or cDNA "
                                        + "coordinate if using -n/--noncoding argument)" )
                                .create( "c" );
      Options options = new Options();
      options.addOption(gene);
      options.addOption(species);
      options.addOption(coordinate);
      options.addOption(b37);
      options.addOption(noncoding);
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
      formatter.printHelp("CodingToGenomic "
              + "[-g <gene>] [-c <coordinate>] [options]\n", o, false);
      System.out.println("");
      System.exit(exitVal);
  }
  
  static class OptionComparator<T extends Option> implements Comparator<T> {

        private static final List<String> OPT_ORDER = Arrays.asList(
                "g", "t", "c", "n", "s", "b", "h");

        @Override
        public int compare(T o1, T o2) {
            return OPT_ORDER.indexOf(o1.getOpt()) - OPT_ORDER.indexOf(o2.getOpt());
        }
    }
}
