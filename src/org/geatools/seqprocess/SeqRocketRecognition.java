package org.geatools.seqprocess;
/*******************************************************************************
 *  ========================================================================
 *  GEATools : a free Genomic Event Analysis Tools for the Java(tm) platform
 *  ========================================================================
 *
 *  (C) Copyright 2016, by Qi Wang and Contributors.
 *
 *  This file is part of GEATools.
 *
 *  GEATools is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 *  [Oracle and Java are registered trademarks of Oracle and/or its affiliates. 
 *  Other names may be trademarks of their respective owners.]
 *
 *  (C) Copyright 2000-2014, by Original Author and Contributors.
 *
 *  Original Author: Qi Wang;
 *  Contributor(s): 
 *  Changes (from 1-Jan-2016)
 *  ---------------------------------------------------------------------------
 *******************************************************************************/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.geatools.data.structure.BlastInfo;
import org.geatools.data.structure.SeqCompo;
import org.geatools.data.structure.SeqCompoAlignInfo;
import org.geatools.data.structure.SeqCompoFeatures;
import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqRocket;
import org.geatools.data.structure.SeqRocketPair;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;

public class SeqRocketRecognition {    
	
	String SEQ_LEFT_SIDE="left";
	String SEQ_RIGHT_SIDE="right";
	
	public static String BARCODE_NAME_DEFINITION="Barcode";
	public static String PRIMER_NAME_DEFINITION="Primer";
	public static String PRIMERCONT_NAME_DEFINITION="PrimerCont";
	public static String BAIT_NAME_DEFINITION="Bait";
	public static String BAITBRK_NAME_DEFINITION="BaitBrk";
	public static String BAITARM_NAME_DEFINITION="BaitArm";
	public static String BODY_NAME_DEFINITION="Body";
	public static String RC3ADAPTER_NAME_DEFINITION="RC3Adapter";
	public static String RC3TOGETHER_NAME_DEFINITION="RC3Together";
	public static String RC3PRIMERCONT_NAME_DEFINITION="RC3PrimerCont";
	public static String RC3PRIMER_NAME_DEFINITION="RC3Primer";
	public static String RC3BARCODE_NAME_DEFINITION="RC3Barcode";
	public static String RC3SEQLINKER_NAME_DEFINITION="RC3SeqLinker";
		  
	public static String BARCODE_COLOR_DEFINITION="black";
	public static String PRIMER_COLOR_DEFINITION="red";
	public static String PRIMERCONT_COLOR_DEFINITION="yellow";
	public static String BAIT_COLOR_DEFINITION="pink";
	public static String BAITBRK_COLOR_DEFINITION="gray";
	public static String BAITARM_COLOR_DEFINITION="brown";
	public static String BODY_COLOR_DEFINITION="green";
	public static String RC3ADAPTER_COLOR_DEFINITION="orange";
	public static String RC3PRIMERCONT_COLOR_DEFINITION="orange";
	public static String RC3PRIMER_COLOR_DEFINITION="blue";
	public static String RC3BARCODE_COLOR_DEFINITION="black";
	public static String RC3SEQLINKER_COLOR_DEFINITION="cyan";
	public static String RC3_COLOR_DEFINITION="blue";

	double territoryLeftExtendRatio=0.3d;
	double territoryRightExtendRatio=1.0d;
	float leftTerritoryPercent=0.6f;
	float rightTerritoryPercent=0.6f;
	
	int barcodeMinLen=8;
	int primerMinLen=8;
	
	boolean doSeqLeft=true;
	boolean doSeqRight=true;
	boolean doSeqRightTogether=true;	 
	
	static String regexSeqID = "(\\s+SeqID=)(\\d+)(#\\s*)";
	  
	static String homeDir=".";
	String dataDir;
	List<String> tmpFiles=new ArrayList<String>();
	String tmpDir;
	
	List<SeqRocket> seqRockets;
	List<SeqRocketPair> seqPairRockets;
	boolean isRocketsOK=false;	
	
	public void setHomeDir(String dir){
	    homeDir=dir;
	}	 
	
	public void setDataDir(String dir){
	    dataDir=dir;
	}
	
	public void setTmpDir(String dir){
	    tmpDir=dir;
	}	 
	  
	public String getTmpDir(){
	    return tmpDir;
	}	 
	
	public void setTmpFiles(List<String> files){
		tmpFiles=files;
	}

	public List<String> getTmpFiles(){
	    return tmpFiles;
	}
	
	public void setDoSeqLeft(boolean isDoIt){ 
	    doSeqLeft=isDoIt;
	}
	
	public void setDoSeqRight(boolean isDoIt){ 
	    doSeqRight=isDoIt;
	}
	
	public void setDoSeqRightTogether(boolean isDoIt){ 
	    doSeqRightTogether=isDoIt;
	}
	
	public boolean isDoSeqLeft(){ 
	    return doSeqLeft;
	}
	
	public boolean isDoSeqRight(){ 
	    return doSeqRight;
	}
	
	public boolean isDoSeqRightTogether(){ 
	    return doSeqRightTogether;
	}

	public void setBarcodeMinLen(int minLen){ 
		if(minLen>0) barcodeMinLen=minLen;
	}
	
	public int getBarcodeMinLen(){ 
		return barcodeMinLen;
	}
	
	public void setPrimerMinLen(int minLen){ 
		if(minLen>0) primerMinLen=minLen;
	}
	
	public int getPrimerMinLen(){ 
		return primerMinLen;
	}

	public boolean isSeqRocketsOK(){
	    return isRocketsOK;
	}
	public void setSeqRockets(List<SeqRocket> rockets){
	    seqRockets=rockets;
	    if(rockets!=null && rockets.size()>0) isRocketsOK=true;
		else isRocketsOK=false;;
	}
	public List<SeqRocket> getSeqRockets(){
	    return seqRockets;
	}  
	public void setSeqPairRockets(List<SeqRocketPair> pairRockets){
	    seqPairRockets=pairRockets;
		if(pairRockets!=null && pairRockets.size()>0) isRocketsOK=true;
		else isRocketsOK=false;;
	}
	public List<SeqRocketPair> getSeqPairRockets(){
	    return seqPairRockets;
	}
	
	public List<String> getRecognizedFileList(){
	    
		if(seqRockets==null) return null;
		
		List<String> seqFiles=new ArrayList<String>();
	    
		for(SeqRocket rocket:seqRockets) {
			seqFiles.add(rocket.recognizedSeqFile);
		}
				
		return seqFiles;
	}

	public List<String> getRecognizedMaskFileList(){
	    
		if(seqRockets==null) return null;
		List<String> seqFiles=new ArrayList<String>();
	    
		for(SeqRocket rocket:seqRockets) {
			seqFiles.add(rocket.recognizedSeqFileMasked);
		}
				
		return seqFiles;
	}
	
	public List<String> getRecognizedTrimFileList(){
	    
		if(seqRockets==null) return null;
		List<String> seqFiles=new ArrayList<String>();
	    
		for(SeqRocket rocket:seqRockets) {
			seqFiles.add(rocket.recognizedSeqFileTrimmed);
		}
				
		return seqFiles;
	}
	
	public List<ArrayList<String>> getRecognizedPairFileList(){
	    
		if(seqPairRockets==null) return null;
		
		List<ArrayList<String>> seqFiles12=new ArrayList<ArrayList<String>>();
		ArrayList<String> seqFiles=new ArrayList<String>();
		ArrayList<String> seqFiles2=new ArrayList<String>();
		for(SeqRocketPair rocket12:seqPairRockets) {
			seqFiles.add(rocket12.forward.recognizedSeqFile);
			seqFiles2.add(rocket12.reverse.recognizedSeqFile);
		}
		
		seqFiles12.add(seqFiles);
		seqFiles12.add(seqFiles2);
		
		seqFiles=null;
		seqFiles2=null;		
		
		return seqFiles12;
	}

	public List<ArrayList<String>> getRecognizedPairMaskFileList(){
	    
		if(seqPairRockets==null) return null;
		
		List<ArrayList<String>> seqFiles12=new ArrayList<ArrayList<String>>();
		ArrayList<String> seqFiles=new ArrayList<String>();
		ArrayList<String> seqFiles2=new ArrayList<String>();
		for(SeqRocketPair rocket12:seqPairRockets) {
			seqFiles.add(rocket12.forward.recognizedSeqFileMasked);
			seqFiles2.add(rocket12.reverse.recognizedSeqFileMasked);
		}
		
		seqFiles12.add(seqFiles);
		seqFiles12.add(seqFiles2);
		
		seqFiles=null;
		seqFiles2=null;
		
		return seqFiles12;
	}
	
	public List<ArrayList<String>> getRecognizedPairTrimFileList(){
	    
		if(seqPairRockets==null) return null;
		
		List<ArrayList<String>> seqFiles12=new ArrayList<ArrayList<String>>();
		ArrayList<String> seqFiles=new ArrayList<String>();
		ArrayList<String> seqFiles2=new ArrayList<String>();
		for(SeqRocketPair rocket12:seqPairRockets) {
			seqFiles.add(rocket12.forward.recognizedSeqFileTrimmed);
			seqFiles2.add(rocket12.reverse.recognizedSeqFileTrimmed);
		}
		
		seqFiles12.add(seqFiles);
		seqFiles12.add(seqFiles2);
		
		seqFiles=null;
		seqFiles2=null;
		
		return seqFiles12;
	}
	
	public List<SeqCompo> getLibExpSeqCompo(String libExpSeqInfoFile){
		    
		    List<SeqCompo> expSeqCompoList = new ArrayList<SeqCompo>();
		    
		    String expName=null;
		    String barcodeName=null;
			String barcodeSeq=null;	
			String primerSeq=null;	
			String primerContSeq=null;	
			String baitSeq=null;
			String baitArmSeq=null;
			String rc3PrimerContSeq=null;	
			String rc3PrimerSeq=null;	
			String rc3BarcodeSeq=null;	
			String rc3SeqLinkerSeq=null;
			
			String freqCutterSeq=null;
			
			boolean isExpNameOK=false;	
			boolean isPrimerOK=false;
			boolean isBarcodeOK=false;
		    int expIdx=0;
		    List<String> expNameList=new ArrayList<String>();
			try{    
			   BufferedReader br;              
			   br = new BufferedReader(new FileReader(libExpSeqInfoFile));
			   String line;
			   String [] itemSplited;
			   String attrName="";
			   String attrValue="";
			   line = br.readLine();		
			   while(true){ 			   
				  if(line == null) break;
				  if(line.trim().indexOf(">")==0){			     
					 line=br.readLine();			
					 if (line == null) break;
					 expName=null;
					 barcodeName=null;
					 barcodeSeq=null;	
					 primerSeq=null;	
					 primerContSeq=null;	
					 baitSeq=null;
					 baitArmSeq=null;
					 rc3PrimerContSeq=null;	
					 rc3PrimerSeq=null;	
					 rc3BarcodeSeq=null;	
					 rc3SeqLinkerSeq=null;
					 isExpNameOK=false;
					 isBarcodeOK=false;
					 isPrimerOK=false;
					 expIdx=expIdx+1;
					 while(line.indexOf(">")<0 && line.indexOf("=")>0){						
						attrName="";
						attrValue="";
						itemSplited=line.split("=");
						if(itemSplited.length>=2){					  
							attrName=itemSplited[0].trim();
							attrValue=itemSplited[1].trim();
						}else if(itemSplited.length==1){
							attrName=itemSplited[0].trim();
							attrValue="";
						}
						if(attrName.equalsIgnoreCase("ExperimentName")){
						    expName=attrValue;
						    if(expName==null || expName.equals("")){
							    System.out.println("!!! Warning: You have empty 'ExperimentName'");
								System.out.println("We automatically use it's index as experiment name.");
								expName="Exp_"+expIdx; 
						    }
						    
						    expName=expName.replaceAll("\\s+","-");
							   
							if(!expNameList.contains(expName)){ 
								expNameList.add(expName);
							}else{
								System.out.println("!!! Warning: You have repeat 'ExperimentName="+expName+"'.");
								System.out.println("We automatically make it unique by adding experiment index.");
								expName=expName+"_"+expIdx;
								expNameList.add(expName);
							}
							isExpNameOK=true;
					    }else if(attrName.equalsIgnoreCase("BarcodeName")){
					    	barcodeName=attrValue;	
							if(barcodeName!=null && !barcodeName.equals("")) barcodeName=barcodeName.replaceAll("\\s+","-");
						}else if(attrName.equalsIgnoreCase("BarcodeSeq")){
							barcodeSeq=attrValue;
							if(barcodeSeq!=null && !barcodeSeq.equals("")) isBarcodeOK=true;
						}else if(attrName.equalsIgnoreCase("PrimerSeq")){
							primerSeq=attrValue;
							if(primerSeq!=null && !primerSeq.equals("")) isPrimerOK=true;
						}else if(attrName.equalsIgnoreCase("PrimerContSeq")){
							primerContSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("BaitTerritorySeq")){
							baitSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("BaitTerritoryArmSeq")){
							baitArmSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("RC3PrimerContSeq")){
							rc3PrimerContSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("RC3PrimerSeq")){
							rc3PrimerSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("RC3BarcodeSeq")){
							rc3BarcodeSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("RC3SeqLinkerSeq")){
							rc3SeqLinkerSeq=attrValue;	
						}else if(attrName.equalsIgnoreCase("FrequentCutterSeq")){
							freqCutterSeq=attrValue;	
						}						
					
						line=br.readLine(); 
						 
						if (line == null) break;					    
					 }
					 
					 if(isBarcodeOK || isPrimerOK){
						 if(!isExpNameOK){
							System.out.println("!!! Warning: You didn't correctly configure 'ExperimentName' for No."+expIdx);
							System.out.println("We automatically use it's index as experiment name.");
							expName="Exp_"+expIdx; 
							expNameList.add(expName);
						 }
						 SeqCompo seqCompo = new SeqCompo();
						 seqCompo.isActive=true;
						 seqCompo.expName=expName;
						 seqCompo.barcodeName=barcodeName;
						 seqCompo.barcodeSeq=barcodeSeq;
						 seqCompo.primerSeq=primerSeq;
						 seqCompo.primerContSeq=primerContSeq;
						 seqCompo.baitTerritorySeq=baitSeq;
						 seqCompo.baitTerritoryArmSeq=baitArmSeq;
						 seqCompo.rc3PrimerContSeq=rc3PrimerContSeq;
						 seqCompo.rc3PrimerSeq=rc3PrimerSeq;
						 seqCompo.rc3BarcodeSeq=rc3BarcodeSeq;
						 seqCompo.rc3SeqLinkerSeq=rc3SeqLinkerSeq;

						 seqCompo.freqCutterSeq=freqCutterSeq;
					
						 expSeqCompoList.add(seqCompo);
						 seqCompo=null;							
					 }else{
						 System.out.println("!!! Warning: You didn't correctly configure 'BarcodeSeq' or 'PrimerSeq' for No."+expIdx);
						 System.out.println("We just skip it."); 
					 }
						
				  }else{
					 line = br.readLine();
				  }				 
				}           		   
				br.close();
				br=null;
				
			 }catch(IOException e){
			        System.out.println(e);
			 }
					
			 return expSeqCompoList;
		 
	}
	  
	List<SeqRocket> buildRockets(String libExpSeqInfoFile, SeqRocketConfig rocketConfig){
		 
			List<SeqCompo> libSeqRocketList=getLibExpSeqCompo(libExpSeqInfoFile);			
			
			//READ_LEFT_SIDE
			SeqCompoRecognizer barcode;
		    SeqCompoRecognizer primer;
		    SeqCompoRecognizer primerCont;
			SeqCompoRecognizer bait;
			SeqCompoRecognizer baitBrk;
			SeqCompoRecognizer baitArm;
			//READ_RIGHT_SIDE
		    SeqCompoRecognizer rc3PrimerCont;
		    SeqCompoRecognizer rc3Primer;
			SeqCompoRecognizer rc3Barcode;
			SeqCompoRecognizer rc3SeqLinker;
			
			SeqCompoRecognizer rc3Together;
			
			//READ_LEFT_SIDE Seq
			String barcodeSeq=null;	
			String primerSeq=null;	
			String primerContSeq=null;	
			String baitSeq=null;
			String baitArmSeq=null;	
			//READ_RIGHT_SIDE Seq
			String rc3PrimerContSeq=null;	
			String rc3PrimerSeq=null;	
			String rc3BarcodeSeq=null;	
			String rc3SeqLinkerSeq=null;
			
			String rc3TogetherSeq=null;
			
			String rocketName="";			
			List<SeqRocket> rockets =new ArrayList<SeqRocket>();
			SeqRocket libSeqRocket;	
			String minLeftRecognizerSeq="";
			String allRightRecognizerSeq="";
			 
			for(int i=0;i<libSeqRocketList.size();i++){
				barcodeSeq=libSeqRocketList.get(i).barcodeSeq;	
			    primerSeq=libSeqRocketList.get(i).primerSeq;	
			    primerContSeq=libSeqRocketList.get(i).primerContSeq;	
			    baitSeq=libSeqRocketList.get(i).baitTerritorySeq;
			    baitArmSeq=libSeqRocketList.get(i).baitTerritoryArmSeq;
			    rc3PrimerContSeq=libSeqRocketList.get(i).rc3PrimerContSeq;	
			    rc3PrimerSeq=libSeqRocketList.get(i).rc3PrimerSeq;	
			    rc3BarcodeSeq=libSeqRocketList.get(i).rc3BarcodeSeq;	
				rc3SeqLinkerSeq=libSeqRocketList.get(i).rc3SeqLinkerSeq;	
				
				minLeftRecognizerSeq=barcodeSeq+primerSeq+primerContSeq;
				allRightRecognizerSeq=rc3PrimerContSeq+rc3PrimerSeq+rc3BarcodeSeq+rc3SeqLinkerSeq;
				 
				rocketName=libSeqRocketList.get(i).expName;	
				if(libSeqRocketList.get(i).barcodeName!=null 
						&& !libSeqRocketList.get(i).barcodeName.trim().equals("")) {
				   rocketName=rocketName+"."+libSeqRocketList.get(i).barcodeName;
				}
				libSeqRocket=new SeqRocket();
				libSeqRocket.rocketName=rocketName;
				libSeqRocket.expName=libSeqRocketList.get(i).expName;
				libSeqRocket.isActive=libSeqRocketList.get(i).isActive;
			
				libSeqRocket.seqRecognizers=new ArrayList<SeqCompoRecognizer>();
				libSeqRocket.seqCompoFeatures=new SeqCompoFeatures();
				libSeqRocket.seqCompoFeatures.compoNames=new ArrayList<String>();	
				libSeqRocket.seqCompoFeatures.compoColors=new ArrayList<String>();
				libSeqRocket.savedCompoAlignedSeqFiles=new ArrayList<String>();		
				
				if(barcodeSeq!=null && barcodeSeq.length()>0 && (primerSeq==null || primerSeq.length()==0)){
				   if(barcodeSeq.length()<barcodeMinLen){
					    System.out.println(
					    	"Warning: The length of barcode for "
					        +libSeqRocketList.get(i).expName+" is less than "+barcodeMinLen+"!"
					    );				    
				   }
				   barcode=new SeqCompoRecognizer();			 
				   libSeqRocket.seqCompoFeatures.compoNames.add(BARCODE_NAME_DEFINITION);		
				   libSeqRocket.seqCompoFeatures.compoColors.add(BARCODE_COLOR_DEFINITION);		
				   barcode.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(BARCODE_NAME_DEFINITION);			
				   barcode.rawSeq=barcodeSeq;			
				   barcode.seqName=rocketName+"."+BARCODE_NAME_DEFINITION;
				   barcode.side=SEQ_LEFT_SIDE; 
				   libSeqRocket.seqRecognizers.add(barcode);		 
				}else if(primerSeq!=null && primerSeq.length()>0){
				   if(primerSeq.length()<primerMinLen){
				        System.out.println(
				    	  "Warning: The length of primer for "
				          +libSeqRocketList.get(i).expName+" is less than "+primerMinLen+"!"
				        );
				        primerMinLen=primerSeq.length();
				   }
				   barcode=new SeqCompoRecognizer();			 
				   libSeqRocket.seqCompoFeatures.compoNames.add(BARCODE_NAME_DEFINITION);			 
				   libSeqRocket.seqCompoFeatures.compoColors.add(BARCODE_COLOR_DEFINITION);			  
				   barcode.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(BARCODE_NAME_DEFINITION);
				   barcode.rawSeq=barcodeSeq;
				   barcode.seqName=rocketName+"."+BARCODE_NAME_DEFINITION;
				   barcode.side=SEQ_LEFT_SIDE; 
				   libSeqRocket.seqRecognizers.add(barcode);
				   
				   primer=new SeqCompoRecognizer();
				   libSeqRocket.seqCompoFeatures.compoNames.add(PRIMER_NAME_DEFINITION);
				   libSeqRocket.seqCompoFeatures.compoColors.add(PRIMER_COLOR_DEFINITION);
				   primer.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(PRIMER_NAME_DEFINITION);			 	
				   primer.rawSeq=primerSeq;			
				   primer.seqName=rocketName+"."+PRIMER_NAME_DEFINITION;		
				   primer.side=SEQ_LEFT_SIDE; 			 
				   libSeqRocket.seqRecognizers.add(primer);  
				}
				
				if(baitSeq!=null && baitSeq.length()>0){
				   primerCont=new SeqCompoRecognizer();		
				   libSeqRocket.seqCompoFeatures.compoNames.add(PRIMERCONT_NAME_DEFINITION);
				   libSeqRocket.seqCompoFeatures.compoColors.add(PRIMERCONT_COLOR_DEFINITION);
				   primerCont.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(PRIMERCONT_NAME_DEFINITION);
				   primerCont.rawSeq=primerContSeq;	
				   primerCont.seqName=rocketName+"."+PRIMERCONT_NAME_DEFINITION;
				   primerCont.side=SEQ_LEFT_SIDE; 
		           libSeqRocket.seqRecognizers.add(primerCont);
		           
		           bait=new SeqCompoRecognizer();		  
		           libSeqRocket.seqCompoFeatures.compoNames.add(BAIT_NAME_DEFINITION);		  
		           libSeqRocket.seqCompoFeatures.compoColors.add(BAIT_COLOR_DEFINITION);		  
				   bait.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(BAIT_NAME_DEFINITION);		  	
				   bait.rawSeq=baitSeq;		  
				   bait.seqName=rocketName+"."+BAIT_NAME_DEFINITION;
				   bait.side=SEQ_LEFT_SIDE; 
				   libSeqRocket.seqRecognizers.add(bait);
				   
				   baitBrk=new SeqCompoRecognizer();		  
				   libSeqRocket.seqCompoFeatures.compoNames.add(BAITBRK_NAME_DEFINITION);		  
				   libSeqRocket.seqCompoFeatures.compoColors.add(BAITBRK_COLOR_DEFINITION);		  
				   baitBrk.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(BAITBRK_NAME_DEFINITION);		  	
				   baitBrk.rawSeq=baitSeq;		  
				   baitBrk.seqName=rocketName+"."+BAITBRK_NAME_DEFINITION;
				   baitBrk.side=SEQ_LEFT_SIDE; 
				   libSeqRocket.seqRecognizers.add(baitBrk);
				   
				   if(baitArmSeq!=null && baitArmSeq.length()>0){
					  baitArm=new SeqCompoRecognizer();		  
					  libSeqRocket.seqCompoFeatures.compoNames.add(BAITARM_NAME_DEFINITION);		  
					  libSeqRocket.seqCompoFeatures.compoColors.add(BAITARM_COLOR_DEFINITION);		  
					  baitArm.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(BAITARM_NAME_DEFINITION);		  	
					  baitArm.rawSeq=baitArmSeq;		  
					  baitArm.seqName=rocketName+"."+BAITARM_NAME_DEFINITION;	
					  baitArm.side=SEQ_LEFT_SIDE; 
					  libSeqRocket.seqRecognizers.add(baitArm);			  
				   }

				}else if(primerContSeq!=null && primerContSeq.length()>0){
				   primerCont=new SeqCompoRecognizer();
				   libSeqRocket.seqCompoFeatures.compoNames.add(PRIMERCONT_NAME_DEFINITION);
				   libSeqRocket.seqCompoFeatures.compoColors.add(PRIMERCONT_COLOR_DEFINITION);
				   primerCont.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(PRIMERCONT_NAME_DEFINITION);	
				   primerCont.rawSeq=primerContSeq;
				   primerCont.seqName=rocketName+"."+PRIMERCONT_NAME_DEFINITION;	
				   primerCont.side=SEQ_LEFT_SIDE; 
				   libSeqRocket.seqRecognizers.add(primerCont);		     		  
				}	
				
				if(doSeqRightTogether) {
				   rc3TogetherSeq="";
				   if(rc3PrimerContSeq!=null && rc3PrimerContSeq.trim().length()>0){
					  rc3TogetherSeq=rc3TogetherSeq+rc3PrimerContSeq;	     		  
				   }	
				
				   if(rc3PrimerSeq!=null && rc3PrimerSeq.trim().length()>0){
					  rc3TogetherSeq=rc3TogetherSeq+rc3PrimerSeq;	  	     		  
				   }
				
				   if(rc3BarcodeSeq!=null && rc3BarcodeSeq.trim().length()>0){
					  rc3TogetherSeq=rc3TogetherSeq+rc3BarcodeSeq;	    		  
				   }
				
				   if(rc3SeqLinkerSeq!=null && rc3SeqLinkerSeq.trim().length()>0){
					  rc3TogetherSeq=rc3TogetherSeq+rc3SeqLinkerSeq;	 
				   }
				   
				   if(rc3TogetherSeq!=null && !rc3TogetherSeq.equals("")){
				     rc3Together=new SeqCompoRecognizer();
				     libSeqRocket.seqCompoFeatures.compoNames.add(RC3TOGETHER_NAME_DEFINITION);
				     libSeqRocket.seqCompoFeatures.compoColors.add(RC3_COLOR_DEFINITION);
				     rc3Together.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(RC3TOGETHER_NAME_DEFINITION);	
				     rc3Together.rawSeq=rc3TogetherSeq;
				     rc3Together.seqName=rocketName+"."+RC3TOGETHER_NAME_DEFINITION;	
				     rc3Together.side=SEQ_RIGHT_SIDE; 
				     libSeqRocket.seqRecognizers.add(rc3Together);
				   }
				}else {
				   if(rc3PrimerContSeq!=null && rc3PrimerContSeq.trim().length()>0){
					  rc3PrimerCont=new SeqCompoRecognizer();
					  libSeqRocket.seqCompoFeatures.compoNames.add(RC3PRIMERCONT_NAME_DEFINITION);
					  libSeqRocket.seqCompoFeatures.compoColors.add(RC3PRIMERCONT_COLOR_DEFINITION);
					  rc3PrimerCont.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(RC3PRIMERCONT_NAME_DEFINITION);	
					  rc3PrimerCont.rawSeq=rc3PrimerContSeq;
					  rc3PrimerCont.seqName=rocketName+"."+RC3PRIMERCONT_NAME_DEFINITION;	
					  rc3PrimerCont.side=SEQ_RIGHT_SIDE;
					  libSeqRocket.seqRecognizers.add(rc3PrimerCont);		     		  
				   }	
						
				   if(rc3PrimerSeq!=null && rc3PrimerSeq.trim().length()>0){
					  rc3Primer=new SeqCompoRecognizer();
					  libSeqRocket.seqCompoFeatures.compoNames.add(RC3PRIMER_NAME_DEFINITION);
					  libSeqRocket.seqCompoFeatures.compoColors.add(RC3PRIMER_COLOR_DEFINITION);
					  rc3Primer.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(RC3PRIMER_NAME_DEFINITION);	
					  rc3Primer.rawSeq=rc3PrimerSeq;
					  rc3Primer.seqName=rocketName+"."+RC3PRIMER_NAME_DEFINITION;	
					  rc3Primer.side=SEQ_RIGHT_SIDE;	  
					  libSeqRocket.seqRecognizers.add(rc3Primer);		     		  
				   }
						
				   if(rc3BarcodeSeq!=null && rc3BarcodeSeq.trim().length()>0){
					  rc3Barcode=new SeqCompoRecognizer();
					  libSeqRocket.seqCompoFeatures.compoNames.add(RC3BARCODE_NAME_DEFINITION);
					  libSeqRocket.seqCompoFeatures.compoColors.add(RC3BARCODE_COLOR_DEFINITION);
					  rc3Barcode.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(RC3BARCODE_NAME_DEFINITION);	
					  rc3Barcode.rawSeq=rc3BarcodeSeq;
					  rc3Barcode.seqName=rocketName+"."+RC3BARCODE_NAME_DEFINITION;	
					  rc3Barcode.side=SEQ_RIGHT_SIDE;		  
					  libSeqRocket.seqRecognizers.add(rc3Barcode);		     		  
				   }
						
				   if(rc3SeqLinkerSeq!=null && rc3SeqLinkerSeq.trim().length()>0){
					  rc3SeqLinker=new SeqCompoRecognizer();
					  libSeqRocket.seqCompoFeatures.compoNames.add(RC3SEQLINKER_NAME_DEFINITION);
					  libSeqRocket.seqCompoFeatures.compoColors.add(RC3SEQLINKER_COLOR_DEFINITION);
					  rc3SeqLinker.index=libSeqRocket.seqCompoFeatures.compoNames.indexOf(RC3SEQLINKER_NAME_DEFINITION);	
					  rc3SeqLinker.rawSeq=rc3SeqLinkerSeq;
					  rc3SeqLinker.seqName=rocketName+"."+RC3SEQLINKER_NAME_DEFINITION;	
					  rc3SeqLinker.side=SEQ_RIGHT_SIDE; 
					  libSeqRocket.seqRecognizers.add(rc3SeqLinker);		     		  
				   }				
				}
				
				libSeqRocket.minLeftRecognizerSeq=minLeftRecognizerSeq;
				libSeqRocket.allRightRecognizerSeq=allRightRecognizerSeq;
				
				rockets.add(libSeqRocket);
				
				libSeqRocket=null;
				
				barcode=null;
				primer=null;
				primerCont=null;
			    bait=null;
				baitBrk=null;
				baitArm=null;		
			
				rc3Together=null;
				
				rc3PrimerCont=null;
				rc3Primer=null;
				rc3Barcode=null;
				rc3SeqLinker=null;
			}
			
			libSeqRocketList=null;
			
			rocketConfig.configRecognizer(rockets);
			
			return rockets;
    }
	  
    
	List<SeqCompoAlignInfo> getSeqObj(String inSeqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=null;
			
		    if(SeqOperation.isFASTASeq(inSeqFile)){		
			   seqObjList=getSeqObjFromFASTA(inSeqFile);	   
			}else if(SeqOperation.isFASTQSeq(inSeqFile)){		
			   seqObjList=getSeqObjFromFASTQ(inSeqFile);	   
			}
		    
			return seqObjList;
	}
		  
	List<SeqCompoAlignInfo> getSeqObjFromFASTA(String seqFile){
			  
			    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
			    SeqCompoAlignInfo perSeq=new SeqCompoAlignInfo();
				int seqID=0;
				try{    
					BufferedReader br;              
			        br = new BufferedReader(new FileReader(seqFile));
				    String line;
					String seqIdentifier;
					String seqLine;
					String seqName;
					String [] itemSplited;
					seqObjList=new ArrayList<SeqCompoAlignInfo>();
					seqID=0;
					line = br.readLine();
					if (line == null){ 
					   br.close();
					   return seqObjList;
					}
					while(line.length()==0 || line.matches("\\s*")){
					   line = br.readLine();
					   if (line == null){ 
						   br.close();
						   return seqObjList;
					   }
					}
					while(true){           
				       if (line == null) break;
				       line=line.trim();
					   if(line.indexOf(">")==0){					    
					     seqIdentifier=line.substring(1,line.length());				 
						 itemSplited=seqIdentifier.split("\\s+");
						 seqName=itemSplited[0].trim();
						 			 
						 line=br.readLine();
						 seqLine="";
						 if (line == null) break;
						 while(line.indexOf(">")<0){
					      seqLine = seqLine+line.trim();
						  line=br.readLine();
						  if (line == null) break;
						 }
						 	  
						 perSeq=new SeqCompoAlignInfo();
						 perSeq.seqID=seqID;
						 perSeq.seqIdentifier=seqIdentifier;
						 perSeq.seqName=seqName;
						 perSeq.seqLength=seqLine.length();
						 perSeq.seq=seqLine;			 
						 seqObjList.add(perSeq);
						 perSeq=null;  
						 
						 seqID=seqID+1;
					   }			 
					}//while
			        		  
					br.close();
					
				}catch(IOException e){
			        System.out.println(e);
			    }
			    //System.out.println("invalid reads :"+(seqNum-uniSeqName.size()));
				//uniSeqName=null;
			    return seqObjList; 
	}
			  
	List<SeqCompoAlignInfo> getSeqObjFromFASTQ(String seqFile){
			  
			    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
			    SeqCompoAlignInfo perSeq=new SeqCompoAlignInfo();
				int seqID=0;
				try{    
					BufferedReader br;              
			        br = new BufferedReader(new FileReader(seqFile));
				    String line;
					String seqIdentifier;	
					String seqLine;
					String seqQualityLine;
					String seqName;
					String [] itemSplited;
					seqObjList=new ArrayList<SeqCompoAlignInfo>();
					seqID=0;
					outerloop:
					while(true){		         	
			           
					   //for read header/identifier
			           line = br.readLine();			
				       if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   line=line.trim();
			           if(line.indexOf("@")!=0) {
					      System.out.println("Error in reading fastq line: "+line);
						  break;
					   }			   
			      	   seqIdentifier=line.substring(1,line.length());				 
					   itemSplited=seqIdentifier.split("\\s+");
					   seqName=itemSplited[0].trim();
					   
					   //for read sequence
					   line=br.readLine();			 
					   if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   seqLine=line.trim();			 
					   //seqLine=seqLine.replaceAll("N","n");
					   
					   //for read '+' tag
					   line=br.readLine();			  
					   if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   line=line.trim();
					   if(line.indexOf("+")!=0) {
						  System.out.println("Error in reading fastq line: "+line);
						  break;
					   }		      
			           
					   //for read base quality encode
			           line=br.readLine();			 
					   if (line == null) break;
					   while(line.length()==0 || line.matches("\\s*")){
						  line = br.readLine();
						  if (line == null) break outerloop;
					   }
					   seqQualityLine=line.trim();
			           
			           perSeq=new SeqCompoAlignInfo();
			           perSeq.seqID=seqID;
			           perSeq.seqIdentifier=seqIdentifier;
					   perSeq.seqName=seqName;
					   perSeq.seqLength=seqLine.length();
					   perSeq.seq=seqLine;				 
				       perSeq.seqQualityEncode=seqQualityLine;			 
					   seqObjList.add(perSeq);
					   perSeq=null;     
					   
					   seqID=seqID+1;
			    
					}//while
			        	  
					br.close();
					
					line=null;
					seqIdentifier=null;
					seqLine=null;
					seqName=null;
					itemSplited=null;
					
				}catch(IOException e){
			        System.out.println(e);
			    }
			    //System.out.println("invalid reads :"+(seqNum-uniSeqName.size()));
				//uniSeqName=null;
			    return seqObjList; 
	}
	
	boolean checkStruct(){
		  
	       boolean isOK=true;
			
		   return isOK; 		   
	}
		  
	List<SeqCompoAlignInfo> checkSeq(String inSeqFile){
			    
			List<SeqCompoAlignInfo> seqObjList = null;

			if(SeqOperation.isFASTASeq(inSeqFile)){			  
			   seqObjList=checkFASTASeq(inSeqFile);			
			}else if(SeqOperation.isFASTQSeq(inSeqFile)){		  
			   seqObjList=checkFASTQSeq(inSeqFile);			
			}
				
			return seqObjList;
	}
		 
	List<SeqCompoAlignInfo> checkFASTASeq(String seqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
		    SeqCompoAlignInfo perSeq;
			int seqID=0;

			try{    

				BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier;
				String seqLine;
				String seqName;
				String [] itemSplited;			
				
				line = br.readLine();
				if (line == null){ 
				   br.close();
				   return seqObjList;
				}
				while(line.length()==0 || line.matches("\\s*")){
				   line = br.readLine();
				   if (line == null){ 
					   br.close();
					   return seqObjList;
				   }
				}
				
				seqID=0;
				while(true){           
			       if (line == null) break;
			       line=line.trim();
				   if(line.indexOf(">")==0){				    
				     seqIdentifier=line.substring(1,line.length());			 
					 itemSplited=seqIdentifier.split("\\s+");
					 seqName=itemSplited[0].trim();
					 seqIdentifier=seqName; ////
					 
					 line=br.readLine();
					 seqLine="";
					 if (line == null) break;
					 while(line.indexOf(">")<0){
				       seqLine = seqLine+line.trim();
					   line=br.readLine();
					   if (line == null) break;
					 }
					 seqLine=seqLine.replaceAll("N","n");					
					 
					 perSeq=new SeqCompoAlignInfo();
					 perSeq.seqID=seqID;
					 perSeq.seqIdentifier=seqIdentifier;
					 perSeq.seqName=seqName;
					 perSeq.seqLength=seqLine.length();
					 perSeq.seq=seqLine;
					 seqObjList.add(perSeq);
					 perSeq=null;   
					 
					 seqID=seqID+1;
				   }
				}           		   
				br.close();
			
				line=null;
				seqIdentifier=null;	
				seqLine=null;
				seqName=null;
				itemSplited=null;
			}catch(IOException e){
		        System.out.println(e);
		    }

		    return seqObjList; 
	}
		 
	List<SeqCompoAlignInfo> checkFASTQSeq(String seqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
		    SeqCompoAlignInfo perSeq;
			int seqID=0;
			
			try{    
		       	BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier;	
				String seqLine;
				String seqName;		
				String seqQualityLine;
				String [] itemSplited;
				seqID=0;
				outerloop:	
				while(true){  		           
				   //for read header/identifier
		           line = br.readLine();			
			       if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
		           if(line.indexOf("@")!=0) {
				      System.out.println("Error in reading fastq line: "+line);
					  break;
				   }			   
		      	   seqIdentifier=line.substring(1,line.length());				 
				   itemSplited=seqIdentifier.split("\\s+");
				   seqName=itemSplited[0].trim();
				   seqIdentifier=seqName; ////
				   
				   //for read sequence
				   line=br.readLine();			 
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   seqLine=line.trim();			 
				   seqLine=seqLine.replaceAll("N","n");
				   
				   //for read '+' tag
				   line=br.readLine();			  
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
				   if(line.indexOf("+")!=0) {
					  System.out.println("Error in reading fastq line: "+line);
					  break;
				   }		      
		           
				   //for read base quality encode
		           line=br.readLine();			 
				   if (line == null) break;
				   while(line.length()==0 || line.matches("\\s*")){
					  line = br.readLine();
					  if (line == null) break outerloop;
				   }
				   seqQualityLine=line.trim();

				   perSeq=new SeqCompoAlignInfo();
				   perSeq.seqID=seqID;
				   perSeq.seqIdentifier=seqIdentifier;
				   perSeq.seqName=seqName;
				   perSeq.seqLength=seqLine.length();
				   perSeq.seq=seqLine;
				   perSeq.seqQualityEncode=seqQualityLine;
				   seqObjList.add(perSeq);
				   perSeq=null;	
				   
				   seqID=seqID+1;	
				}           		   
				br.close();
			
		    	line=null;
		    	seqIdentifier=null;	
				seqLine=null;
				seqName=null;
				itemSplited=null;
			}catch(IOException e){
		        System.out.println(e);
		    }

		    return seqObjList; 
	}
	
	void createLeftRecognizerSeq(List<SeqCompoAlignInfo> seqObjList, SeqRocket seqRocket, 
			  SeqCompoRecognizer recognizer, String outSeqFile){
	  
		try{    
	        int trimSEndIndex=-1;
			int trimLeftShift=0;
			int trimSEnd=0;
			for(int i=0;i<seqRocket.seqRecognizers.size();i++){
			  if(seqRocket.seqRecognizers.get(i).done 
					  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_LEFT_SIDE) 
					  && seqRocket.seqRecognizers.get(i).leftShiftForNext){		 
			     trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
			     trimLeftShift=seqRocket.seqRecognizers.get(i).trimLeftShift;
			  }
			}
			
	        int territoryLen=recognizer.territoryLen;
				
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));
			String nameLine;
			String seqLine="";
			String trimmedSeqLine="";
		    int seqID=0;
			int seqLen=0;
	    	if(trimSEndIndex>=0){
			  if(recognizer.leftSubForBLAST){
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;	
				  trimmedSeqLine="nnnnnnnnnn";			
				  if(trimSEnd+territoryLen<seqLen)
					 trimmedSeqLine=seqLine.substring(trimSEnd,trimSEnd+territoryLen);
				  else if(trimSEnd<seqLen)
					 trimmedSeqLine=seqLine.substring(trimSEnd,seqLen);
					 
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;	
				  trimmedSeqLine="nnnnnnnnnn";	
				  if(trimSEnd<seqLen) trimmedSeqLine=seqLine.substring(trimSEnd,seqLen);			 
				  
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		 
			  }
			}else{
			  if(recognizer.leftSubForBLAST){
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
					
				  if(territoryLen<seqLen)
					trimmedSeqLine=seqLine.substring(0,territoryLen);
				  else 
					trimmedSeqLine=seqLine.substring(0,seqLen);
					 
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();			 
				  
				  writer.write(seqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		  
			  }
			
			}
			writer.close();
			
		}catch(IOException e){
	        System.out.println(e);
	    }  
	 
	}
	
	void createRightRecognizerSeq(List<SeqCompoAlignInfo> seqObjList, SeqRocket seqRocket, 
			  SeqCompoRecognizer recognizer, String outSeqFile){
	  
		try{    
	        int trimSEndIndex=-1;
			int trimRightShift=0;
			int trimSEnd=0;
			for(int i=seqRocket.seqRecognizers.size()-1;i>=0;i--){
			  if(seqRocket.seqRecognizers.get(i).done 
					  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_RIGHT_SIDE)
					  && seqRocket.seqRecognizers.get(i).rightShiftForNext){		 
			     trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
			     trimRightShift=seqRocket.seqRecognizers.get(i).trimRightShift;
			  }
			}
			
	        int territoryLen;
				
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));
			String nameLine;
			String seqLine="";
			String trimmedSeqLine="";
		    int seqID=0;
			int seqLen=0;
	    	if(trimSEndIndex>=0){
			  if(recognizer.rightSubForBLAST){
				territoryLen=recognizer.territoryLen;
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)+trimRightShift;
				  if(trimSEnd>seqLen) trimSEnd=seqLen;
				  trimmedSeqLine="nnnnnnnnnn";			
				  if(trimSEnd-territoryLen>=0)
					trimmedSeqLine=seqLine.substring(trimSEnd-territoryLen,trimSEnd);
				  else 
					trimmedSeqLine=seqLine.substring(0,trimSEnd);
					 
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)+trimRightShift;
				  if(trimSEnd>seqLen) trimSEnd=seqLen;
				  trimmedSeqLine="nnnnnnnnnn";	
				  if(trimSEnd<seqLen) trimmedSeqLine=seqLine.substring(0,trimSEnd);			 
				  
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		 
			  }
			}else{
			  if(recognizer.rightSubForBLAST){				
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  territoryLen=Math.max(recognizer.territoryLen,Math.round(seqLen*recognizer.territoryPercent));	
				  if(territoryLen<seqLen) {
					trimmedSeqLine=seqLine.substring(seqLen-territoryLen,seqLen);
				  }else { 
					trimmedSeqLine=seqLine.substring(0,seqLen);
				  }	 
				  writer.write(trimmedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();			 
				  
				  writer.write(seqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		  
			  }
			
			}
			writer.close();
			
		}catch(IOException e){
	        System.out.println(e);
	    }  
	 
	}

	
	void createBLASTTarSeq(List<SeqCompoAlignInfo> seqObj,String outSeqFile){
  		
		try{    

			String seqLine;
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
		      
				 seqLine=seqObj.get(i).seq;
				 writer.write(">"+i);
				 writer.newLine();
				 writer.write(seqLine);
				 writer.newLine();
				 writer.flush(); 
				 
			}           		   
			
			writer.close();
		}catch(IOException e){
	        System.out.println(e);
	    }

    } 

	  
	void createLeftSubBLASTTarSeq(List<SeqCompoAlignInfo> seqObj,int leftSubSeqLen,String outSeqFile){
		  		
			try{    

				String seqLine;
				//String seqName;
				int subSeqLength=0;
				BufferedWriter writer=null;
		        writer=new BufferedWriter(new FileWriter(outSeqFile));	
				for(int i=0;i<seqObj.size();i++){     
			      
					 seqLine=seqObj.get(i).seq;
					 subSeqLength=leftSubSeqLen;
					 if(leftSubSeqLen>seqLine.length()) subSeqLength=seqLine.length();
		             seqLine=seqLine.substring(0,subSeqLength);
					  
					 writer.write(">"+i);
					 writer.newLine();
					 writer.write(seqLine);
					 writer.newLine();
					 writer.flush(); 
					 
				}           		   
				
				writer.close();
			}catch(IOException e){
		        System.out.println(e);
		    }

	} 
	
	void createRightSubBLASTTarSeq(List<SeqCompoAlignInfo> seqObj,int rightSubSeqLen,String outSeqFile){
  		
		try{    

			String seqLine;
			//String seqName;
			int subSeqLength=0;
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
		      
				 seqLine=seqObj.get(i).seq;
				 subSeqLength=rightSubSeqLen;
				 if(rightSubSeqLen>seqLine.length()) subSeqLength=seqLine.length();
	             seqLine=seqLine.substring(seqLine.length()-subSeqLength,seqLine.length());
				  
				 writer.write(">"+i);
				 writer.newLine();
				 writer.write(seqLine);
				 writer.newLine();
				 writer.flush(); 
				 
			}           		   
			
			writer.close();
		}catch(IOException e){
	        System.out.println(e);
	    }

    } 


	void initSeqAlignArray(List<SeqCompoAlignInfo> seqObjList, int seqCompoNum){  			

		      ArrayList<Integer> seqAlignSStart;	
		      ArrayList<Integer> seqAlignSEnd;  
		    
		      for(int s=0;s<seqObjList.size();s++){   
			        
					 seqAlignSStart=new ArrayList<Integer>();
					 seqAlignSEnd=new ArrayList<Integer>();
			         for(int i=0;i<seqCompoNum;i++){	           
						seqAlignSStart.add(-1);
						seqAlignSEnd.add(-1);
			         }
			         seqObjList.get(s).seqAlignSStart=seqAlignSStart;
			         seqObjList.get(s).seqAlignSEnd=seqAlignSEnd;
		             seqAlignSStart=null;
					 seqAlignSEnd=null;		
				   		 
			  }	//for	
		     
	}
	
	List<SeqCompoAlignInfo> getRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, int compoIndx){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;		
		int alignSStartIndex=compoIndx;		  
		int alignSEndIndex=alignSStartIndex;	
		
		for(int i=0; i<seqObjList.size(); i++){
		
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)>=0 
				   && seqObjList.get(i).seqAlignSEnd.get(alignSEndIndex)>=0){
			   
		     seqObj.add(seqObjList.get(i));
		     
		   }
		}
		
		return seqObj; 
	}
	 
	List<SeqCompoAlignInfo> getNoRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, int compoIndx){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;	
		
		int alignSStartIndex=compoIndx;		  
		//int alignSEndIndex=alignSStartIndex;
		
		for(int i=0; i<seqObjList.size(); i++){
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)==-1)
		     seqObj.add(seqObjList.get(i));
		
		}
		
		return seqObj;
	 
	} 


	List<SeqCompoAlignInfo> getRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, 
			SeqRocket seqRocket,String seqCompoName){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;		
		int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);		  
		int alignSEndIndex=alignSStartIndex;	
		
		for(int i=0; i<seqObjList.size(); i++){
		
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)>=0 
				   && seqObjList.get(i).seqAlignSEnd.get(alignSEndIndex)>=0){
			   
		     seqObj.add(seqObjList.get(i));
		     
		   }
		}
		
		return seqObj; 
	}
	 
	List<SeqCompoAlignInfo> getNoRecognizedSeq(List<SeqCompoAlignInfo> seqObjList,
			SeqRocket seqRocket,String seqCompoName){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;	
		
		int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);		  
		//int alignSEndIndex=alignSStartIndex;
		
		for(int i=0; i<seqObjList.size(); i++){
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)==-1)
		     seqObj.add(seqObjList.get(i));		
		}
		
		return seqObj;
	 
	} 
	  
	void setSeqLeftExactAlignInfo(List<SeqCompoAlignInfo> seqObjList,SeqRocket seqRocket,
			String seqCompoName){
	  
	 	try{
			String seqLine;		
			int leftSStartIndex=0;
	        int leftSStart=0;
	        int leftSEnd=0;		
			int trimSEndIndex=-1;
			int trimLeftShift=0;
			int trimSEnd=0;		
			for(int i=0;i<seqRocket.seqRecognizers.size();i++){
			  if(seqRocket.seqRecognizers.get(i).done 
					  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_LEFT_SIDE) 
					  && seqRocket.seqRecognizers.get(i).leftShiftForNext){		 
				  trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
			      trimLeftShift=seqRocket.seqRecognizers.get(i).trimLeftShift;
			  }
			}
			
			int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);		  
			int alignSEndIndex=alignSStartIndex;
			
			SeqCompoRecognizer recognizer=seqRocket.seqRecognizers.get(
					 seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName));	
			String recogniSeq=recognizer.seq;
	        int leftMaxStartIndex=recognizer.exactMaxStart-1;
	        if (leftMaxStartIndex<0) leftMaxStartIndex=0;
	        
	        if(trimSEndIndex>=0){
			  for(int i=0;i<seqObjList.size();i++){
				seqLine=seqObjList.get(i).seq;				
				leftSStartIndex=seqLine.indexOf(recogniSeq);
			    trimSEnd=seqObjList.get(i).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   	
				if(leftSStartIndex>=trimSEnd && leftSStartIndex<=trimSEnd+leftMaxStartIndex){
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			  }
		   }else{
			  for(int i=0;i<seqObjList.size();i++){
				seqLine=seqObjList.get(i).seq;					
				leftSStartIndex=seqLine.indexOf(recogniSeq);		   	
				if(leftSStartIndex>=0 && leftSStartIndex<=leftMaxStartIndex){
				//if(redSStartIndex>=0){	 
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			  }
		   }		
	    }catch(Exception e){
	        System.out.println(e);
	    }

	}
	
	
	void setSeqRightExactAlignInfo(List<SeqCompoAlignInfo> seqObjList,SeqRocket seqRocket,
			String seqCompoName){
		  
	 	try{
			String seqLine;		
			int leftSStartIndex=0;
	        int leftSStart=0;
	        int leftSEnd=0;		
			int trimSStartIndex=-1;
			int trimRightShift=0;
			int trimSStart=0;		
			for(int i=seqRocket.seqRecognizers.size()-1;i>=0;i--){
			  if(seqRocket.seqRecognizers.get(i).done 
					  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_RIGHT_SIDE) 
					  && seqRocket.seqRecognizers.get(i).rightShiftForNext){		 
				  trimSStartIndex=seqRocket.seqRecognizers.get(i).index;
				  trimRightShift=seqRocket.seqRecognizers.get(i).trimRightShift;
			  }
			}
			
			int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);		  
			int alignSEndIndex=alignSStartIndex;
			
			SeqCompoRecognizer recognizer=seqRocket.seqRecognizers.get(
					 seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName));	
			String recogniSeq=recognizer.seq;
	        
			int maxStart;
			int minStart;
	        if(trimSStartIndex>=0){
			  for(int i=0;i<seqObjList.size();i++){
				seqLine=seqObjList.get(i).seq;				
				leftSStartIndex=seqLine.indexOf(recogniSeq);
			    trimSStart=seqObjList.get(i).seqAlignSStart.get(trimSStartIndex)+trimRightShift;	
			    maxStart=trimSStart-recogniSeq.length();
			    minStart=Math.max(0,trimSStart-recognizer.territoryLen);
				if(leftSStartIndex<=maxStart && leftSStartIndex>=minStart){
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			  }
		   }else{			 
			  for(int i=0;i<seqObjList.size();i++){
				seqLine=seqObjList.get(i).seq;					
				leftSStartIndex=seqLine.indexOf(recogniSeq);
				maxStart=seqLine.length()-recogniSeq.length();
				minStart=Math.min(
					Math.round(seqLine.length()*(1-recognizer.territoryPercent)),
					seqLine.length()-recognizer.territoryLen
				);
				minStart=Math.max(0,minStart);
				if(leftSStartIndex<=maxStart && leftSStartIndex>=minStart){				 
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			  }
		   }		
	    }catch(Exception e){
	        System.out.println(e);
	    }

	}

    List<SeqCompoAlignInfo> getLeftSideBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
   		 SeqCompoRecognizer recognizer, String seqCompoName, SeqRocket seqRocket, String blastOutFile){
   	
   		List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> ();
   		SeqCompoAlignInfo seq;
   		
   		int alignLen=0;
   		int mismatchNum=0;
   		int gapNum=0;
   		int qStart=0;
   		int sStart=0;
   		int sEnd=0;
   		//int sEnd0=0;
   		
   		int minAlignLen=recognizer.minAlignLen;
   		int maxMismatchNum=recognizer.maxMismatchNum;
   		int maxGapNum=recognizer.maxGapNum;
   		int maxQStart=recognizer.maxQStart;
   		int maxSStart=recognizer.maxSStart;	
   		double maxMismatchRatio=recognizer.maxMismatchRatio;
   		double maxGapRatio=recognizer.maxGapRatio;
   		
   		//ArrayList<ArrayList <String>> blastOutList=new ArrayList<ArrayList <String>> ();
   	    List<ArrayList <String>> blastOut=FileOperation.getMatrixFromFile(blastOutFile);
   		List <Integer> seqID=new ArrayList <Integer>();
   		List <Integer> sStartList=new ArrayList <Integer>();
   		List <Integer> sEndList=new ArrayList <Integer>();
   		String readName="";
   		int readID0=-1;
   		int readID=0;
   		
   		for(int i=0; i<blastOut.size();i++){
   		  readName=blastOut.get(i).get(BlastInfo.colSName).trim();
   		  readID=Integer.parseInt(readName);
   		  if(readID>readID0){
   				
   			 alignLen=Integer.parseInt(blastOut.get(i).get(BlastInfo.colAlignLen).trim());
   			 mismatchNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colMismatchNum).trim());
   			 gapNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colGapNum).trim());
   			 qStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQStart).trim());
   			 sStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSStart).trim());
   			 sEnd=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSEnd).trim());	
   	         /*
   			 if(sStart>sEnd){
   			   sEnd0=sStart;
   			   sStart=sEnd;
   			   sEnd=sEnd0;
   			 }
   			 */
   	         maxMismatchNum=(int) Math.round(alignLen*maxMismatchRatio);
   			 maxGapNum=(int) Math.round(alignLen*maxGapRatio);	 			
   				
   			 if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum && gapNum<=maxGapNum 
   						&& qStart<=maxQStart && sStart<=maxSStart && sEnd>=sStart){
   				seqID.add(readID); 
   				sStartList.add(sStart);
   				sEndList.add(sEnd);
   				readID0=readID;
   			    //blastOutList.add(blastOut.get(i));
   			}
   		  }
   		}
   	
   		blastOut=null;		
   		//blastOutList=null;
   		  
   		int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);	
   		int alignSEndIndex=alignSStartIndex;	    
   		
   		int trimSEndIndex=-1;
   		int trimLeftShift=0;
   		int trimSEnd=0;
   		for(int i=0;i<seqRocket.seqRecognizers.size();i++){
   		  if(seqRocket.seqRecognizers.get(i).done 
   				  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_LEFT_SIDE)
   				  && seqRocket.seqRecognizers.get(i).leftShiftForNext ){		 
   			
   			  trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
   			  trimLeftShift=seqRocket.seqRecognizers.get(i).trimLeftShift;
   		  
   		  }
   		}
   			
   		if(trimSEndIndex>=0){
   			for(int i=0;i<seqID.size();i++){
   			   seq =new SeqCompoAlignInfo();	
   			   seq=seqObjList.get(seqID.get(i));
   			   trimSEnd=seq.seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   		   
   			   seq.seqAlignSStart.set(alignSStartIndex,trimSEnd+sStartList.get(i));
   			   seq.seqAlignSEnd.set(alignSEndIndex,trimSEnd+sEndList.get(i));
   			   seqObj.add(seq);
   			   seq=null;
   			}
   		}else{
   			for(int i=0;i<seqID.size();i++){
   			   seq =new SeqCompoAlignInfo();	
   			   seq=seqObjList.get(seqID.get(i));
   			   seq.seqAlignSStart.set(alignSStartIndex,sStartList.get(i));
   			   seq.seqAlignSEnd.set(alignSEndIndex,sEndList.get(i));
   			   seqObj.add(seq);
   			   seq=null;
   			}	
   		}

   		seqID=null;
   		sStartList=null;
   		sEndList=null;
   		
   		return seqObj;
   	  
    } 

    List<ArrayList<SeqCompoAlignInfo>> getBaitBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
   		 SeqCompoRecognizer recognizer, String seqCompoName, SeqRocket seqRocket, String blastOutFile){
   	
   		List<ArrayList<SeqCompoAlignInfo>>  bait=new  ArrayList<ArrayList<SeqCompoAlignInfo>> ();
   		ArrayList<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> ();
   		ArrayList<SeqCompoAlignInfo> noSeqObj=new ArrayList<SeqCompoAlignInfo> ();
   		SeqCompoAlignInfo seq;
   		
   		int alignLen=0;
   		int mismatchNum=0;
   		int gapNum=0;
   		int qStart=0;
   		int qEnd=0;
   		int sStart=0;
   		int sEnd=0;
   		int sEnd0=0;
   		
   		int minAlignLen=recognizer.minAlignLen;
   		int maxMismatchNum=recognizer.maxMismatchNum;
   		int maxGapNum=recognizer.maxGapNum;
   		int maxQStart=recognizer.maxQStart;
   		int maxSStart=recognizer.maxSStart;	
   		double maxMismatchRatio=recognizer.maxMismatchRatio;
   		double maxGapRatio=recognizer.maxGapRatio;
   		int baitBrkQPos=0;
   			
   	    List<ArrayList <String>> blastOut=FileOperation.getMatrixFromFile(blastOutFile);   	
   		List <Integer> sStartList=new ArrayList <Integer>();
   		List <Integer> sEndList=new ArrayList <Integer>();
   		List <Integer> baitBrkQPosList=new ArrayList <Integer>();
   		List <Integer> BLASTSeqIDList=new ArrayList <Integer>();
   		List <Integer> noBLASTSeqIDList=new ArrayList <Integer>();
   		String readName="";
   		int readID0=-1;
   		int readID=0;
   		
   		for(int i=0; i<blastOut.size();i++){
   		  readName=blastOut.get(i).get(BlastInfo.colSName).trim();
   		  readID=Integer.parseInt(readName);
   		  if(readID>readID0){			
   			 alignLen=Integer.parseInt(blastOut.get(i).get(BlastInfo.colAlignLen).trim());
   			 mismatchNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colMismatchNum).trim());
   			 gapNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colGapNum).trim());
   			 qStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQStart).trim());
   			 qEnd=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQEnd).trim());
   			 sStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSStart).trim());
   			 sEnd=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSEnd).trim());	
   	         if(sStart>sEnd){
   			    sEnd0=sStart;
   				sStart=sEnd;
   				sEnd=sEnd0;
   			 }
   	         maxMismatchNum=(int) Math.round(alignLen*maxMismatchRatio);
   			 maxGapNum=(int) Math.round(alignLen*maxGapRatio);
   			 if(seqCompoName.equals(BAIT_NAME_DEFINITION)){
   				if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum 
   							&& gapNum<=maxGapNum && qStart<=maxQStart && sStart<=maxSStart){
   				   BLASTSeqIDList.add(readID); 
   				   sStartList.add(sStart);
   				   sEndList.add(sEnd);				  
   				   baitBrkQPosList.add(qEnd);
   				   readID0=readID;				
   				}
   			 }else if(seqCompoName.equals(BAITBRK_NAME_DEFINITION)){
   				baitBrkQPos=seqObjList.get(readID).baitBrkQPos;	
   				maxSStart=recognizer.maxSStart-baitBrkQPos;
   				if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum && gapNum<=maxGapNum 
   						&& qStart>baitBrkQPos && qStart<=maxQStart && sStart<=maxSStart){
   				   BLASTSeqIDList.add(readID); 
   				   sStartList.add(sStart);
   				   sEndList.add(sEnd);
   				   readID0=readID;					
   				}   				
   			 }else if(seqCompoName.equals(BAITARM_NAME_DEFINITION)){
   				if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum 
   							&& gapNum<=maxGapNum && qStart<=maxQStart && sStart<=maxSStart){
   				   BLASTSeqIDList.add(readID); 
   				   sStartList.add(sStart);
   				   sEndList.add(sEnd);				  
   				   readID0=readID;				
   				}
   			 }
   		  }
   		}
   	
   		blastOut=null;		
   		  
   		int alignSStartIndex=0;	
   		alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);	
   		int alignSEndIndex=alignSStartIndex;	    
   		
   		int trimSEndIndex=-1;
   		int trimLeftShift=0;
   		int trimSEnd=0;
   		for(int i=0;i<seqRocket.seqRecognizers.size();i++){
   		  if(seqRocket.seqRecognizers.get(i).done 
   				  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_LEFT_SIDE)
   				  && seqRocket.seqRecognizers.get(i).leftShiftForNext ){		 
   		     trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
   		     trimLeftShift=seqRocket.seqRecognizers.get(i).trimLeftShift;
   		  }
   		}
   		
   		if(trimSEndIndex>=0){
   		  for(int i=0;i<BLASTSeqIDList.size();i++){
   			 seq =new SeqCompoAlignInfo();	
   			 seq=seqObjList.get(BLASTSeqIDList.get(i));
   			 trimSEnd=seq.seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   		   
   			 seq.seqAlignSStart.set(alignSStartIndex,trimSEnd+sStartList.get(i));
   			 seq.seqAlignSEnd.set(alignSEndIndex,trimSEnd+sEndList.get(i));
   			 if(seqCompoName.equals(BAIT_NAME_DEFINITION)){
   			    seq.baitBrkQPos=baitBrkQPosList.get(i);
   			    seq.baitBrkSPos=trimSEnd+sEndList.get(i);
   			 }
   			 seqObj.add(seq);
   			 seq=null;
   		  }
   		}else{
   		  for(int i=0;i<BLASTSeqIDList.size();i++){
   			 seq =new SeqCompoAlignInfo();	
   			 seq=seqObjList.get(BLASTSeqIDList.get(i));
   			 seq.seqAlignSStart.set(alignSStartIndex,sStartList.get(i));
   			 seq.seqAlignSEnd.set(alignSEndIndex,sEndList.get(i));
   			 if(seqCompoName.equals(BAIT_NAME_DEFINITION)){
   			    seq.baitBrkQPos=baitBrkQPosList.get(i);
   			    seq.baitBrkSPos=sEndList.get(i);
   			 }
   			 seqObj.add(seq);
   			 seq=null;
   		  }	
   		}
   		
   		bait.add(seqObj);
   		seqObj=null;
   		
   	    //for seq without BLAST alignment
   		for(int i=0; i<seqObjList.size();i++){
   		  if(!BLASTSeqIDList.contains(i)) {
   			noBLASTSeqIDList.add(i);
   			//System.out.println(i);
   		  }
   		}
   		BLASTSeqIDList=null;
   		sStartList=null;
   		sEndList=null;
   		baitBrkQPosList=null;   		
   		for(int i=0;i<noBLASTSeqIDList.size();i++){
   		    seq =new SeqCompoAlignInfo();	
   			seq=seqObjList.get(noBLASTSeqIDList.get(i));
   			noSeqObj.add(seq);
   			seq=null;
   		}   		
   		bait.add(noSeqObj);   		
   		noSeqObj=null;
   		noBLASTSeqIDList=null;
   		
   		return bait;
   	  
    }
    
    List<ArrayList<SeqCompoAlignInfo>> getRightSideBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
      		 SeqCompoRecognizer recognizer, String seqCompoName, SeqRocket seqRocket, String blastOutFile){
         		
	   	List<ArrayList<SeqCompoAlignInfo>> rightSide=new  ArrayList<ArrayList<SeqCompoAlignInfo>> ();
	   	ArrayList<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> ();
	   	ArrayList<SeqCompoAlignInfo> noSeqObj=new ArrayList<SeqCompoAlignInfo> ();
	   	SeqCompoAlignInfo seq;
      		
      	int alignLen=0;
      	int mismatchNum=0;
      	int gapNum=0;
      	int qStart=0;
      	int sStart=0;
      	int sEnd=0;
      
      	int readLen=0;            		
        int territoryLen=recognizer.territoryLen;    
      	int maxMismatchNum=recognizer.maxMismatchNum;
      	int maxGapNum=recognizer.maxGapNum;
      	int minAlignLen=recognizer.minAlignLen;      	
      	double maxMismatchRatio=recognizer.maxMismatchRatio;
      	double maxGapRatio=recognizer.maxGapRatio;   
      	double minAlignRatio=recognizer.minAlignRatio;
      	int maxQStart=0;
      	int minQStart=0;
      	int maxSStart=0;	
      	int minSStart=0;
 
      	List<ArrayList <String>> blastOut=FileOperation.getMatrixFromFile(blastOutFile);      	
      	List <Integer> sStartList=new ArrayList <Integer>();
      	List <Integer> sEndList=new ArrayList <Integer>();
      	List <Integer> BLASTSeqIDList=new ArrayList <Integer>();
       	List <Integer> noBLASTSeqIDList=new ArrayList <Integer>();
      		
      	String readName="";
      	int readID0=-1;
      	int readID=0;
      		
      	for(int i=0; i<blastOut.size();i++){
      		readName=blastOut.get(i).get(BlastInfo.colSName).trim();
      		readID=Integer.parseInt(readName);
      		
      		if(readID>readID0){
      				
      		   alignLen=Integer.parseInt(blastOut.get(i).get(BlastInfo.colAlignLen).trim());
      		   mismatchNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colMismatchNum).trim());
      		   gapNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colGapNum).trim());
      		   qStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQStart).trim());
      		   sStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSStart).trim());
      		   sEnd=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSEnd).trim());	
      	       
      		   readLen=seqObjList.get(readID).seq.length();       		       		   
   			 
      		   if(recognizer.seqLength<=(readLen-sStart+1)){ 
      		     minAlignLen=Math.max(minAlignLen,(int) (recognizer.seqLength*minAlignRatio));
      		   }
      		   if(recognizer.seqLength<=minAlignLen) minAlignLen=recognizer.seqLength;  
 			   
      		   minSStart=1;
      		   if(!recognizer.rightSubForBLAST) { 
      			 minSStart=(int) (0.6f*seqRocket.minLeftRecognizerSeq.length()+1);
      		   }
      		   maxSStart=readLen-alignLen+1;      		       
      		   maxQStart=recognizer.maxQStart;
      		   maxQStart=Math.min(maxQStart, recognizer.seqLength-minAlignLen+1);
      		   maxMismatchNum=(int) Math.round(alignLen*maxMismatchRatio);
      		   maxGapNum=(int) Math.round(alignLen*maxGapRatio);      			

      		   if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum && gapNum<=maxGapNum 
      					&& qStart>=minQStart && qStart<=maxQStart 
      					&& sStart>=minSStart && sStart<=maxSStart && sEnd>=sStart){
      				BLASTSeqIDList.add(readID); 
      				sStartList.add(sStart);
      				sEndList.add(sEnd);
      				readID0=readID;
      			    //blastOutList.add(blastOut.get(i));
      		   }  	
      		}

      	}
      	
      	blastOut=null;		
      	//blastOutList=null;
      		  
      	int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompoName);	
      	int alignSEndIndex=alignSStartIndex;	    
      		
      	int trimSEndIndex=-1;
      	int trimLSideRightShift=0;
      	int trimSEnd=0;
      	int sBaseline=0;
      	for(int i=seqRocket.seqRecognizers.size()-1;i>=0;i--){
      	   if(seqRocket.seqRecognizers.get(i).done 
      				  && seqRocket.seqRecognizers.get(i).side.equalsIgnoreCase(SEQ_RIGHT_SIDE)
      				  && seqRocket.seqRecognizers.get(i).rightShiftForNext){		 
      			
      		  trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
      		  trimLSideRightShift=seqRocket.seqRecognizers.get(i).trimRightShift;
      		  
      	   }
      	}
      			
      	if(trimSEndIndex>=0){
      	   for(int i=0;i<BLASTSeqIDList.size();i++){
      		   seq =new SeqCompoAlignInfo();	
      		   seq=seqObjList.get(BLASTSeqIDList.get(i));
      		   trimSEnd=seq.seqAlignSEnd.get(trimSEndIndex)+trimLSideRightShift;  
      		   sBaseline=0;
      		   if(recognizer.rightSubForBLAST) sBaseline=Math.max(0,trimSEnd-territoryLen);
      		   seq.seqAlignSStart.set(alignSStartIndex,sBaseline+sStartList.get(i));
      		   seq.seqAlignSEnd.set(alignSEndIndex,sBaseline+sEndList.get(i));
      		   seqObj.add(seq);
      		   seq=null;
      	   }
      	}else{
      	   for(int i=0;i<BLASTSeqIDList.size();i++){
      		   seq =new SeqCompoAlignInfo();	
      		   seq=seqObjList.get(BLASTSeqIDList.get(i));
      		   sBaseline=0;
      		   if(recognizer.rightSubForBLAST) sBaseline=Math.max(0,seq.seq.length()-territoryLen);
      		   seq.seqAlignSStart.set(alignSStartIndex,sBaseline+sStartList.get(i));
      		   seq.seqAlignSEnd.set(alignSEndIndex,sBaseline+sEndList.get(i));
      		   seqObj.add(seq);
      		   seq=null;
      	   }	
      	}      		
       	rightSide.add(seqObj);
       	seqObj=null;
       		
       		//for seq without BLAST alignment
       	for(int i=0; i<seqObjList.size();i++){
       	   if(!BLASTSeqIDList.contains(i)) {
       		  noBLASTSeqIDList.add(i);
       			//System.out.println(i);
       	   }
       	}
       	BLASTSeqIDList=null;
       	sStartList=null;
       	sEndList=null;       		
       	for(int i=0;i<noBLASTSeqIDList.size();i++){
       	   seq =new SeqCompoAlignInfo();	
       	   seq=seqObjList.get(noBLASTSeqIDList.get(i));
       	   noSeqObj.add(seq);
       	   seq=null;
       	}       		
       	rightSide.add(noSeqObj);       		
       	noSeqObj=null;
       	noBLASTSeqIDList=null;
      		
      	return rightSide;
      	  
    } 
    
    
    static boolean checkPairEndSeq(List<SeqCompoAlignInfo>seqObjList, List<SeqCompoAlignInfo>seqObjList2){
		
		 boolean isOK=true;
		 if(seqObjList.size()!=seqObjList2.size()){
			System.out.println("Error: different seq num of pair-end seq");
			return false;
		 }
	
		 for(int i=0;i<seqObjList.size();i++){
			if(!seqObjList.get(i).seqName.trim().equals(seqObjList2.get(i).seqName.trim())){ 
				isOK=false;
				System.out.println("Error: inconsistent seq name between pair-end seq");
				System.out.println("Forward: "+seqObjList.get(i).seqName
						+"\t Reverse:"+seqObjList2.get(i).seqName);
				break;
			}
		 }	
		 
		 return isOK;
    }
 
	int getSeqCompoMaxTerritoryLen(List<SeqRocket> rockets, String seqCompoName){		
			
		int len=0;
		int maxTerritoryLen=0;
		for(SeqRocket rocket:rockets) {
			len=rocket.seqRecognizers.get(
					rocket.seqCompoFeatures.compoNames.indexOf(seqCompoName)
			  ).territoryLen;
				
			if(len>maxTerritoryLen) maxTerritoryLen=len;
		}
			
		return maxTerritoryLen;
	 }	
	
	
	 List<SeqCompoAlignInfo> getNonExactAndSaveExactSeq(List<SeqCompoAlignInfo> seqObj,
			 List<SeqRocket> rockets, String seqCompoName, String inSeqFileFormat){
		
		 List<SeqCompoAlignInfo> nonExactMatchedSeqList=new ArrayList<SeqCompoAlignInfo>();
		 try{    
		
			boolean saveAsFASTAFormat=false;
			boolean saveAsFASTQFormat=false;
			
			if(inSeqFileFormat.equalsIgnoreCase("fasta") 
					|| inSeqFileFormat.equalsIgnoreCase("fna") 
					|| inSeqFileFormat.equalsIgnoreCase("fa")){
			   saveAsFASTAFormat=true;
			}else if(inSeqFileFormat.equalsIgnoreCase("fastq") 
					|| inSeqFileFormat.equalsIgnoreCase("fnq") 
					|| inSeqFileFormat.equalsIgnoreCase("fq")){
			   saveAsFASTQFormat=true;
			}
			
			String seqLine;
			String seqIdentifier;	
			String seqQuality;
			String recognizerSeq="";
			boolean isMatch=false;	
		    int matchStartIndex=0;
			String exactAlignedSeqFile="";
				
			SeqCompoRecognizer recognizer=null;
			ArrayList<BufferedWriter> writerList=new ArrayList<BufferedWriter> ();
			BufferedWriter writer0=null;
			for(int j=0;j<rockets.size();j++){
			  exactAlignedSeqFile=rockets.get(j).seqRecognizers.get(
				         rockets.get(j).seqCompoFeatures.compoNames.indexOf(seqCompoName)
				      ).exactAlignedSeqFile;
			  rockets.get(j).seqRecognizers.get(
					     rockets.get(j).seqCompoFeatures.compoNames.indexOf(seqCompoName)
					  ).exactAlignedSeqFile=exactAlignedSeqFile+"."+inSeqFileFormat;
			  writer0=new BufferedWriter(new FileWriter(exactAlignedSeqFile+"."+inSeqFileFormat));	
			  writerList.add(writer0);
		      writer0=null;	
			}
		
			if(saveAsFASTAFormat){
				for(int i=0;i<seqObj.size();i++){     
				   seqIdentifier=seqObj.get(i).seqIdentifier;
				   seqLine=seqObj.get(i).seq;			
				   isMatch=false;
				   for(int j=0;j<rockets.size();j++){
					  recognizer=rockets.get(j).seqRecognizers.get(
							 rockets.get(j).seqCompoFeatures.compoNames.indexOf(seqCompoName));
					  recognizerSeq=recognizer.seq;
					  matchStartIndex=seqLine.indexOf(recognizerSeq);			 
					  if(matchStartIndex>=0 && matchStartIndex<=recognizer.exactMaxStart-1){
						isMatch=true;				
						writerList.get(j).write(">"+seqIdentifier);
						writerList.get(j).newLine();
						writerList.get(j).write(seqLine);
						writerList.get(j).newLine();
						writerList.get(j).flush(); 
					  } 
				   }
					
				   if(!isMatch){					 
				      nonExactMatchedSeqList.add(seqObj.get(i));		
				   }
				}
			}else if(saveAsFASTQFormat){
				for(int i=0;i<seqObj.size();i++){     
				   seqIdentifier=seqObj.get(i).seqIdentifier;
				   seqLine=seqObj.get(i).seq;			
				   seqQuality=seqObj.get(i).seqQualityEncode;
					
				   isMatch=false;
				   for(int j=0;j<rockets.size();j++){
					  recognizer=rockets.get(j).seqRecognizers.get(
							    rockets.get(j).seqCompoFeatures.compoNames.indexOf(seqCompoName)
							 );
					  recognizerSeq=recognizer.seq;
					  matchStartIndex=seqLine.indexOf(recognizerSeq);			 
					  if(matchStartIndex>=0 && matchStartIndex<=recognizer.exactMaxStart-1){
						isMatch=true;			
						writerList.get(j).write("@"+seqIdentifier);
						writerList.get(j).newLine();
						writerList.get(j).write(seqLine);
						writerList.get(j).newLine();
						writerList.get(j).write("+");
						writerList.get(j).newLine();
						writerList.get(j).write(seqQuality);
						writerList.get(j).newLine();
						writerList.get(j).flush(); 
					  } 
				   }
				   if(!isMatch){					 
					  nonExactMatchedSeqList.add(seqObj.get(i));		
				   }
				}
			}
			
		    for(int j=0;j<writerList.size();j++){		
			    writerList.get(j).close();
			}
			//writer.close();
			
			seqLine=null;
			seqIdentifier=null;		
			recognizerSeq=null;
			
		 }catch(IOException e){
		    System.out.println(e);
		 }
		
		 return nonExactMatchedSeqList;
	
	 } 
	 
	 List<String> getCompoExactAlignedFiles(List<SeqRocket> rockets,String seqCompoName){
		  
	       List<String> exactFiles=new ArrayList<String>();
			
			for(int f=0;f<rockets.size();f++){
				exactFiles.add(
				  rockets.get(f).seqRecognizers.get(
					    rockets.get(f).seqCompoFeatures.compoNames.indexOf(seqCompoName)
				  ).exactAlignedSeqFile
				);
			}
			
			return exactFiles;
	 }
	 
	 void saveRecognizedSeq(List<SeqCompoAlignInfo> seqObj,SeqRocket seqRocket){ 	    
	
			int trimREndLeftShift=0;
			int trimLEndRightShift=0;	
		
			int leftTrimSEndIndex=-1;
			int leftMaskSEndIndex=-1;
			int rightTrimSStartIndex=-1;
			int rightMaskSStartIndex=-1;		
		
	        boolean leftTrimSave=false;
			boolean leftMaskSave=false;
			boolean rightTrimSave=false;
			boolean rightMaskSave=false;
			
			SeqCompoRecognizer recognizer;
			SeqCompoRecognizer leftRecognizer=null;
			SeqCompoRecognizer leftTrimRecognizer=null;
			SeqCompoRecognizer leftMaskRecognizer=null;
			
			for(int i=0;i<seqRocket.seqRecognizers.size();i++){
			  recognizer=seqRocket.seqRecognizers.get(i);		 
			  if(recognizer.done && recognizer.side.equalsIgnoreCase(SEQ_LEFT_SIDE)){
				leftRecognizer=recognizer; 
			    if(recognizer.leftTrimSave){	
				  leftTrimSave=true;
				  leftTrimSEndIndex=recognizer.index;	
				  trimREndLeftShift=recognizer.trimLeftShift;
				  leftTrimRecognizer=recognizer;
			    }
			    
			    if(recognizer.leftMaskSave){	
				  leftMaskSave=true;
				  leftMaskSEndIndex=recognizer.index;	
				  leftMaskRecognizer=recognizer;
			    }
			  }		  
			}
			
			for(int i=seqRocket.seqRecognizers.size()-1;i>=0;i--){
			  recognizer=seqRocket.seqRecognizers.get(i);
			  if(recognizer.done && recognizer.side.equalsIgnoreCase(SEQ_RIGHT_SIDE)){				  
			    if(recognizer.rightTrimSave){	
				  rightTrimSave=true;
				  rightTrimSStartIndex=recognizer.index;	
				  trimLEndRightShift=recognizer.trimRightShift;		
				}
				    
				if(recognizer.rightMaskSave){	
				  rightMaskSave=true;
				  rightMaskSStartIndex=recognizer.index;				
				}
			  }		  
			}
			
	       if(leftRecognizer==null) return;
	       
			try{
				
			  if(seqRocket.saveRecognizedSeqAsFASTQ){
				 if(leftRecognizer.isFASTQSaved){
					seqRocket.recognizedSeqFile=leftRecognizer.alignedSeqFile;
				 }else{
				    leftRecognizer.alignedSeqFile=leftRecognizer.alignedSeqFile+".fastq";
				    saveAsFASTQFile(seqObj,leftRecognizer.alignedSeqFile);
				    seqRocket.recognizedSeqFile=leftRecognizer.alignedSeqFile;
				 }
				 seqRocket.isRecognized=true;
			  }else if(seqRocket.saveRecognizedSeqAsFASTA){
				 if(leftRecognizer.isFASTASaved){
					seqRocket.recognizedSeqFile=leftRecognizer.alignedSeqFile;
				 }else{
				    leftRecognizer.alignedSeqFile=leftRecognizer.alignedSeqFile+".fna";
				    saveAsFASTAFile(seqObj,leftRecognizer.alignedSeqFile);
				    seqRocket.recognizedSeqFile=leftRecognizer.alignedSeqFile;
				 }
				 seqRocket.isRecognized=true;
			  }
				
			  String seqLine;
			  String seqIdentifier;			
			  String trimmedSeq="";       
		      int trimIndex=0;
		      BufferedWriter writer=null; 
		                
		      //// save trimmed seq
			  seqRocket.recognizedSeqFileTrimmed=null;
			  seqRocket.isEndTrimmed=false;		 
			  if(leftTrimSave && rightTrimSave){
				  writer=new BufferedWriter(new FileWriter(leftTrimRecognizer.alignedSeqFileTrimmed));	
				  for(int i=0;i<seqObj.size();i++){     
					seqIdentifier=seqObj.get(i).seqIdentifier;
					seqLine=seqObj.get(i).seq;				
					//trim right	
					trimIndex=seqObj.get(i).seqAlignSStart.get(rightTrimSStartIndex)+trimLEndRightShift;	
					if(trimIndex>=0 && trimIndex<seqLine.length()) {				
					   trimmedSeq=seqLine.substring(0,trimIndex);
				       if (trimmedSeq.length()==0) trimmedSeq="nnnnnnnnnn";
					
					   seqLine=trimmedSeq;
					}
					
					//trim left
					trimIndex=seqObj.get(i).seqAlignSEnd.get(leftTrimSEndIndex)-trimREndLeftShift;
					if(trimIndex<0) trimIndex=0;
					if(trimIndex<seqLine.length()) 
					  trimmedSeq=seqLine.substring(trimIndex,seqLine.length());
					else
					  trimmedSeq="nnnnnnnnnn";
						
					if (trimmedSeq.length()==0) trimmedSeq="nnnnnnnnnn";
						
					writer.write(">"+seqIdentifier);
					writer.newLine();
					writer.write(trimmedSeq);
					writer.newLine();
					writer.flush();                 				
				  }
				  writer.close();
				  seqRocket.recognizedSeqFileTrimmed=leftTrimRecognizer.alignedSeqFileTrimmed;
				  seqRocket.isEndTrimmed=true;
			  }else if(leftTrimSave && !rightTrimSave){
				  writer=new BufferedWriter(new FileWriter(leftTrimRecognizer.alignedSeqFileTrimmed));	
				  for(int i=0;i<seqObj.size();i++){     
					 seqIdentifier=seqObj.get(i).seqIdentifier;
					 seqLine=seqObj.get(i).seq;	
						  
					 trimIndex=seqObj.get(i).seqAlignSEnd.get(leftTrimSEndIndex)-trimREndLeftShift;
					 if(trimIndex<0) trimIndex=0;			
					 if(trimIndex<seqLine.length()) 
					   trimmedSeq=seqLine.substring(trimIndex,seqLine.length());
					 else
					   trimmedSeq="nnnnnnnnnnnnn";
						
					 if (trimmedSeq.length()==0) trimmedSeq="nnnnnnnnnn";
						
					 writer.write(">"+seqIdentifier);
					 writer.newLine();
					 writer.write(trimmedSeq);
					 writer.newLine();
					 writer.flush();                 				
				  }
				  writer.close();
				  seqRocket.recognizedSeqFileTrimmed=leftTrimRecognizer.alignedSeqFileTrimmed;
				  seqRocket.isEndTrimmed=true;
			  }else if(rightTrimSave && !leftTrimSave){
				  writer=new BufferedWriter(new FileWriter(leftRecognizer.alignedSeqFileTrimmed));	
				  for(int i=0;i<seqObj.size();i++){     
					 seqIdentifier=seqObj.get(i).seqIdentifier;
					 seqLine=seqObj.get(i).seq;			
					 trimmedSeq=seqLine;  
					 trimIndex=seqObj.get(i).seqAlignSStart.get(rightTrimSStartIndex)+trimLEndRightShift;				
					 if(trimIndex>=0 && trimIndex<seqLine.length()) {				
					    trimmedSeq=seqLine.substring(0,trimIndex);
					 }
						
					 if (trimmedSeq.length()==0) trimmedSeq="nnnnnnnnnn";
						
					 writer.write(">"+seqIdentifier);
					 writer.newLine();
					 writer.write(trimmedSeq);
					 writer.newLine();
					 writer.flush();                				
				  }
				  writer.close();
				  seqRocket.recognizedSeqFileTrimmed=leftRecognizer.alignedSeqFileTrimmed;
				  seqRocket.isEndTrimmed=true;
			  } 
	         
			  //// save masked seq
			  String maskedSeq="";
			  int maskPos=0;
			  StringBuilder str;
			  String tmp="";
			  seqRocket.recognizedSeqFileMasked=null;
			  seqRocket.isEndMasked=false;
			  if(leftMaskSave && rightMaskSave){
				  writer=new BufferedWriter(new FileWriter(leftMaskRecognizer.alignedSeqFileMasked));	
				  for(int i=0;i<seqObj.size();i++){     
					 seqIdentifier=seqObj.get(i).seqIdentifier;
					 seqLine=seqObj.get(i).seq;	
					 // mask right side
					 maskPos=seqObj.get(i).seqAlignSStart.get(rightMaskSStartIndex);	
		             if(maskPos<0 || maskPos>seqLine.length()) maskPos=seqLine.length()+1;				
									
					 str=new StringBuilder();	
					 for(int x = maskPos-1; x < seqLine.length(); x++){
		                str.append("n");
		             }
					 tmp=str.toString();
					 str=null;				
		             str=new StringBuilder(seqLine);				
					 maskedSeq=str.replace(maskPos-1,seqLine.length(),tmp).toString();
					 seqLine=maskedSeq;
						
					 //mask left side
					 maskPos=seqObj.get(i).seqAlignSEnd.get(leftMaskSEndIndex);
					 if(maskPos<0) maskPos=0;
						
					 str=new StringBuilder();	
					 for (int x = 0; x < maskPos; x++){
		                str.append("n");
		             }
					 tmp=str.toString();
					 str=null;
		             str=new StringBuilder(seqLine);				
					 maskedSeq=str.replace(0,maskPos,tmp).toString();
		             str=null;		
					 tmp=null;
					 writer.write(">"+seqIdentifier);
					 writer.newLine();
					 writer.write(maskedSeq);
					 writer.newLine();
					 writer.flush();                				
				 }
				 writer.close();
				 seqRocket.recognizedSeqFileMasked=leftMaskRecognizer.alignedSeqFileMasked;
				 seqRocket.isEndMasked=true;
			  }else if(leftMaskSave && !rightMaskSave){
				 writer=new BufferedWriter(new FileWriter(leftMaskRecognizer.alignedSeqFileMasked));	
				 for(int i=0;i<seqObj.size();i++){     
				     seqIdentifier=seqObj.get(i).seqIdentifier;
					 seqLine=seqObj.get(i).seq;						  
					 maskPos=seqObj.get(i).seqAlignSEnd.get(leftMaskSEndIndex);
					 if(maskPos<0) maskPos=0;
						
					 str=new StringBuilder();	
					 for (int x = 0; x < maskPos; x++){
		                str.append("n");
		             }
					 tmp=str.toString();
					 str=null;
		             str=new StringBuilder(seqLine);				
					 maskedSeq=str.replace(0,maskPos,tmp).toString();
		             str=null;		
					 tmp=null;
					 writer.write(">"+seqIdentifier);
					 writer.newLine();
					 writer.write(maskedSeq);
					 writer.newLine();
					 writer.flush();                				
				 }
				 writer.close();
				 seqRocket.recognizedSeqFileMasked=leftMaskRecognizer.alignedSeqFileMasked;
				 seqRocket.isEndMasked=true;
			  }else if(rightMaskSave && !leftMaskSave){
				 writer=new BufferedWriter(new FileWriter(leftRecognizer.alignedSeqFileMasked));	
				 for(int i=0;i<seqObj.size();i++){     
					seqIdentifier=seqObj.get(i).seqIdentifier;
					seqLine=seqObj.get(i).seq;				
					  
					maskPos=seqObj.get(i).seqAlignSStart.get(rightMaskSStartIndex);	
		            if(maskPos<0 || maskPos>seqLine.length()) maskPos=seqLine.length()+1;				
									
					str=new StringBuilder();	
					for (int x = maskPos-1; x < seqLine.length(); x++){
		               str.append("n");
		            }
					tmp=str.toString();
					str=null;				
		            str=new StringBuilder(seqLine);				
					maskedSeq=str.replace(maskPos-1,seqLine.length(),tmp).toString();
		            str=null;
					tmp=null;
					writer.write(">"+seqIdentifier);
					writer.newLine();
					writer.write(maskedSeq);
					writer.newLine();
					writer.flush();                				
				 }
				 writer.close();
				 seqRocket.recognizedSeqFileMasked=leftRecognizer.alignedSeqFileMasked;
				 seqRocket.isEndMasked=true;
			  }	
				
			}catch(IOException e){
		      System.out.println(e);
		    }	
	
	 }
	 
	 
	 void saveRecognizedSeqAsHTML(List<SeqCompoAlignInfo> seqObj,SeqRocket seqRocket){ 				
            
		    if(seqObj==null || seqRocket==null) return;
			SeqCompoRecognizer recognizer;
			SeqCompoRecognizer leftRecognizer=null;

			for(int i=0;i<seqRocket.seqRecognizers.size();i++){
			  recognizer=seqRocket.seqRecognizers.get(i);		 
			  if(recognizer.done && recognizer.side.equalsIgnoreCase(SEQ_LEFT_SIDE)){
				leftRecognizer=recognizer; 			    
			  }		  
			}
			
	       if(leftRecognizer==null) return;
	       
		   try{
				
			  String seqLine;
			  String seqIdentifier;			

		      BufferedWriter writer=null; 
				
			  if(seqRocket.saveRecognizedSeqAsHTML){
				 seqRocket.recognizedSeqHTMLFile=leftRecognizer.alignedSeqHTMLFile;
			     String line="";	
			     String line1="";
			     String line2="";
				 int sStart=0;
				 int sEnd=0;
				//int sStart0=0;
				 int sEnd0=0;	
			     //int c=0;			
				 writer=new BufferedWriter(new FileWriter(seqRocket.recognizedSeqHTMLFile));
				 writer.write("<html><title></title><body>");
				 writer.newLine();
				 writer.flush(); 
				 for(int i=0;i<seqObj.size();i++){     
				    seqIdentifier=seqObj.get(i).seqIdentifier;
					seqLine=seqObj.get(i).seq;					 
					
			        writer.write(">"+seqIdentifier+"<br>");
					writer.newLine();
					line="";
					line1="";
					line2="";
					sStart=0;
					sEnd=0;
					//sStart0=0;
					sEnd0=0;			
					for(int k=0;k<seqObj.get(i).seqAlignSStart.size();k++){
					   sStart=seqObj.get(i).seqAlignSStart.get(k);
					   sEnd=seqObj.get(i).seqAlignSEnd.get(k);
					   if(sStart>=0 && sEnd>=0){
					     if(sStart-1>sEnd0){ 
			                line1=seqLine.substring(sEnd0,sStart-1)+
			                      "<span style='color:"+seqRocket.seqCompoFeatures.compoColors.get(k)+"'>";
						 }else{
							sStart=sEnd0+1;
							line1="<span style='color:"+seqRocket.seqCompoFeatures.compoColors.get(k)+"'>";
						 }
							  
						 if(sStart-1<sEnd){
						    line2=seqLine.substring(sStart-1,sEnd)+"</span>";
						 }else{
							   //line2="</span>";
						    if(sStart-1>sEnd0) 
			                  line1=seqLine.substring(sEnd0,sStart-1);
						    else
							  line1="";
							line2="";
						 }
							   
						 line=line+line1+line2;
						 sEnd0=sEnd;				  
					   }				
				    }
					if(sEnd0<seqLine.length()) line=line+seqLine.substring(sEnd0,seqLine.length());
					line=line+"<br>";
						 
					writer.write(line);
					writer.newLine();
					writer.flush(); 				 
				 }
			     writer.write("</body></html>");
				 writer.newLine();
				 writer.flush(); 
						
			     writer.close();
			     
			  }		
				
			}catch(IOException e){
		      System.out.println(e);
		    }	
	
	 }

	
	 void saveAsFASTAFile(List<SeqCompoAlignInfo> seqObj, String outFASTAFile){    	
	
			try{
				String seqLine;
				String seqIdentifier;			
				
				BufferedWriter writer=null;	    
		        writer=new BufferedWriter(new FileWriter(outFASTAFile));
			
				for(int i=0;i<seqObj.size();i++){    
					
				     seqIdentifier=seqObj.get(i).seqIdentifier;
					 seqLine=seqObj.get(i).seq;				 
				
					 writer.write(">"+seqIdentifier);
					 writer.newLine();
					 writer.write(seqLine);
					 writer.newLine();
					 writer.flush(); 		 
				}		
				
		        writer.close();
		        writer=null;	  
			}catch(IOException e){
		        System.out.println(e);
		    }	
	
	 }
	
	
	 void saveAsFASTQFile(List<SeqCompoAlignInfo> seqObj, String outFASTQFile){    	
	
		try{
			String seqLine;
			String seqIdentifier;
	        String seqQuality="";		
			
			BufferedWriter writer=null;	    
	        writer=new BufferedWriter(new FileWriter(outFASTQFile));
		
			for(int i=0;i<seqObj.size();i++){     
			    seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;				 
			    seqQuality=seqObj.get(i).seqQualityEncode;	
				 
	       	    writer.write("@"+seqIdentifier);
				writer.newLine();
				writer.write(seqLine);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush(); 			 
			}		
			
	        writer.close();
	        writer=null;
		}catch(IOException e){
	        System.out.println(e);
	    }	
	
	 }
	
	  
	 void copyAlignInfo(List<SeqCompoAlignInfo> seqObjList, 
			  SeqRocket seqRocket,String fromSeqCompo, String toSeqCompo){
		
		if(fromSeqCompo!=null && toSeqCompo!=null ) { 
		    int fromSStartIndex=0;		  
			fromSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(fromSeqCompo);
			int fromSEndIndex=fromSStartIndex;
		    
		    int toSStartIndex=0;		  
			toSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(toSeqCompo);
			int toSEndIndex=toSStartIndex;		
		
			for(int i=0;i<seqObjList.size();i++){
		
				if(seqObjList.get(i).seqAlignSStart.get(toSStartIndex)==-1 
						&& seqObjList.get(i).seqAlignSEnd.get(toSEndIndex)==-1){
				    
					seqObjList.get(i).seqAlignSStart.set(toSStartIndex,
							seqObjList.get(i).seqAlignSStart.get(fromSStartIndex));
				    seqObjList.get(i).seqAlignSEnd.set(toSEndIndex,
				    		seqObjList.get(i).seqAlignSEnd.get(fromSEndIndex));
				
				}
			}	
		}
		  
	} 	 
	
	List<SeqCompoAlignInfo> combineSeqObj(List<SeqCompoAlignInfo> seqObjList1, 
			 List<SeqCompoAlignInfo> seqObjList2){
	 
	    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo> ();
		
	    if(seqObjList1!=null && seqObjList1.size()>0) seqObjList.addAll(seqObjList1);

	    if(seqObjList2!=null && seqObjList2.size()>0) seqObjList.addAll(seqObjList2);
	    
		return seqObjList;
	}
	
	public static void setSeqNumEncode(List<SeqCompoAlignInfo> seqInfoList, int encodeBaseLen){
	    for(SeqCompoAlignInfo seqInfo:seqInfoList){
	    	setSeqNumEncode(seqInfo,encodeBaseLen);
		}	
	}
	 
	public static void setSeqNumEncode(List<SeqCompoAlignInfo> seqInfoList, float encodeBaseRate){
	    if(encodeBaseRate>1) encodeBaseRate=1;
	    if(encodeBaseRate<0) encodeBaseRate=0;
		int encodeBaseLen;
	    for(SeqCompoAlignInfo seqInfo:seqInfoList){
	       encodeBaseLen=(int) (seqInfo.seq.length()*encodeBaseRate+0.5);
	       setSeqNumEncode(seqInfo,encodeBaseLen);
	    }	
    }

	public static void setSeqNumEncode(SeqCompoAlignInfo seqInfo, int encodeBaseLen){
		 if(encodeBaseLen<1) encodeBaseLen=1;
		 Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int step=seqInfo.seq.length()/encodeBaseLen;	
		 if(step==0) step=1;
		 int i=1;
		 for(int b = 0; b < seqInfo.seq.length(); b=b+step){
			if(i>encodeBaseLen) break;
			i++;
			ch = seqInfo.seq.charAt(b);			   
			val =SeqOperation.getBase2NumMap().get(ch);
			if(val!=null) encodeValStr=encodeValStr+val;		 
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 seqInfo.seqNumEncode=encodeVal;	
	}
	
	/*
	public static void setSeqEncodeVal(SeqCompoAlignInfo seqInfo,int encodeBaseLen){
			 
		 Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int baseLen=Math.min(seqInfo.seq.length(), encodeBaseLen);
		 for (int b = 0; b < baseLen; b++) {
			 ch = seqInfo.seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null) encodeValStr=encodeValStr+val;			 
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 seqInfo.seqNumEncode=encodeVal;
			
	}
	*/		 
	
	public static void setSeqNumRevEncode(List<SeqCompoAlignInfo> seqInfoList, int encodeBaseLen){
	     for(SeqCompoAlignInfo seqInfo:seqInfoList){
	    	 setSeqNumRevEncode(seqInfo,encodeBaseLen);
	     }	
	}
	    
	public static void setSeqNumRevEncode(List<SeqCompoAlignInfo> seqInfoList, float encodeBaseRate){
	     if(encodeBaseRate>1) encodeBaseRate=1;
	     if(encodeBaseRate<0) encodeBaseRate=0;
		 int encodeBaseLen;
	   	 for(SeqCompoAlignInfo seqInfo:seqInfoList){
	    	 encodeBaseLen=(int) (seqInfo.seq.length()*encodeBaseRate+0.5);
	    	 setSeqNumRevEncode(seqInfo,encodeBaseLen);
	     }	
	}

	public static void setSeqNumRevEncode(SeqCompoAlignInfo seqInfo,int encodeBaseLen){
		 if(encodeBaseLen<1) encodeBaseLen=1;
	     Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int seqLen=seqInfo.seq.length();
		 int step=seqInfo.seq.length()/encodeBaseLen;
		 if(step==0) step=1;
		 int i=1;
		 for (int b = seqLen-1; b >=0; b=b-step) {
			 if(i>encodeBaseLen) break;
			 i++;
			 ch = seqInfo.seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null) encodeValStr=encodeValStr+val;			 
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 seqInfo.seqNumRevEncode=encodeVal;		
   }

	 /*
   public static void setSeqNumRevEncode(SeqCompoAlignInfo seqInfo,int encodeBaseLen){
		 
	     Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int seqLen=seqInfo.seq.length();
		 int baseLen=Math.min(seqLen, encodeBaseLen);
		 for (int b = seqLen-1; b > seqLen-baseLen-1; b--) {
			 ch = seqInfo.seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null) encodeValStr=encodeValStr+val;			 
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 seqInfo.seqNumRevEncode=encodeVal;
		
   }
	*/

   
   public static long getSeqNumEncode(String seq, float encodeBaseRate){		 
	     if(encodeBaseRate>1) encodeBaseRate=1;
	     if(encodeBaseRate<0) encodeBaseRate=0;
	     int encodeBaseLen=(int) (seq.length()*encodeBaseRate+0.5);
		 long encodeVal=getSeqNumEncode(seq,encodeBaseLen);		
		 return encodeVal;		
   }
   
   public static long getSeqNumEncode(String seq, int encodeBaseLen){
		 if(encodeBaseLen<1) encodeBaseLen=1;
		 Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int step=seq.length()/encodeBaseLen;
		 if(step==0) step=1;
		 int i=1;
		 for(int b = 0; b < seq.length(); b=b+step){
			if(i>encodeBaseLen) break;
			i++;
			ch = seq.charAt(b);			   
			val =SeqOperation.getBase2NumMap().get(ch);
			if(val!=null) encodeValStr=encodeValStr+val;		 
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 
		 return encodeVal;
   }

  /*
	public static long getSeqNumEncode(String seq,int encodeBaseLen){
		 
	     Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int baseLen=Math.min(seq.length(), encodeBaseLen);
		 for (int b = 0; b < baseLen; b++) {
			 ch = seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null){
				encodeValStr=encodeValStr+val;
			 }
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		
		 return encodeVal;
		
    }
*/
   
	public static long getSeqNumRevEncode(String seq,float encodeBaseRate){		 
	     if(encodeBaseRate>1) encodeBaseRate=1;
	     if(encodeBaseRate<0) encodeBaseRate=0;
		 int encodeBaseLen=(int) (seq.length()*encodeBaseRate+0.5);
		 long encodeVal=getSeqNumRevEncode(seq,encodeBaseLen);		
		 return encodeVal;		
    }
	
	public static long getSeqNumRevEncode(String seq,int encodeBaseLen){		 
		 if(encodeBaseLen<1) encodeBaseLen=1;
		 Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";		 
		 int seqLen=seq.length();
		 int step=seq.length()/encodeBaseLen;	
		 if(step==0) step=1;
		 int i=1;
		 for (int b = seqLen-1; b >=0; b=b-step) {
			 if(i>encodeBaseLen) break;
			 i++;
			 ch = seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null)	encodeValStr=encodeValStr+val;
		 }
		 encodeVal=Long.parseLong(encodeValStr);		 
		
		 return encodeVal;		
   }

	/*
	public static long getSeqNumEncodeRev(String seq,int encodeBaseLen){
		 
	     Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 
		 int seqLen=seq.length();
		 int baseLen=Math.min(seqLen, encodeBaseLen);
		 for (int b = seqLen-1; b > seqLen-baseLen-1; b--) {
			 ch = seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null){
				encodeValStr=encodeValStr+val;
			 }
		 }
		 encodeVal=Long.parseLong(encodeValStr);		 
		
		 return encodeVal;
		
   }
*/
	 
}
