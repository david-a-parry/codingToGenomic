/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.github.codingtogenomic;

/**
 *
 * @author david
 */
public class IdParser {
    private final String id;
    private final boolean isTranscript;
    private final boolean isEnsemblId;
    private boolean isEnsp = false;
    private final String identifierType; 
    IdParser(String s){
        id = s;
        if (id.matches("ENS\\w*G\\d{11}(.\\d+)*")){
            isTranscript = false;
            isEnsemblId = true;
            identifierType = "Ensembl Gene ID";
        }else if (id.matches("ENS\\w*T\\d{11}(.\\d+)*")){
            isTranscript = true;
            isEnsemblId = true;
            identifierType = "Ensembl Transcript ID";
        }else if (id.matches("CCDS\\d+(.\\d+)*")){
            isTranscript = true;
            isEnsemblId = false;
            identifierType = "CCDS ID";
        }else if (id.matches("ENS\\w*P\\d{11}(.\\d+)*")){
            isEnsemblId = true;
            isTranscript = false;
            identifierType = "Ensembl Protein ID";
            isEnsp = true;
        }else if (id.matches("[XN][MPR]_\\d+(.\\d+)*")){
            isTranscript = true;
            isEnsemblId = false;
            identifierType = "RefSeq ID";
        }else if (id.matches("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")){
            isTranscript = true;
            isEnsemblId = false;
            identifierType = "Uniprot ID";
        }else if (id.matches("\\d+")){
            isTranscript = false;
            isEnsemblId = false;
            identifierType = "Entrez Gene ID";
        }else{
            isTranscript = false;
            isEnsemblId = false;
            identifierType = "Gene Symbol or other";
        }
    }
    
    public String getId(){
        return id;
    }
    
    public boolean isTranscript(){
        return isTranscript;
    }
    
    public boolean isEnsp(){
        return isEnsp;
    }
    
    public boolean isEnsemblId(){
        return isEnsemblId;
    }
    
    public String getIdentifierType(){
        return identifierType;
    }
}
