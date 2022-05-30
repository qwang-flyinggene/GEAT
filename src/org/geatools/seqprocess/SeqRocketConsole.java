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

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList; 
import java.util.Calendar;
import java.util.List;

import org.geatools.GEAT;
import org.geatools.data.structure.SeqCompoAlignInfo;
import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqRocket;
import org.geatools.data.structure.SeqRocketPair;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;

public class SeqRocketConsole extends SeqRocketRecognition{   
  
  final public static String LIBSEQINFO_PEMODE_PW="pw"; //pairwise
  final public static String LIBSEQINFO_PEMODE_O2O="o2o"; //one to one
  
  int barcodeMaxTerritoryLen=12;  
  
  public SeqRocketConsole(){ 

  }
  public SeqRocketConsole(String homeDir,String dataDir,String tmpDir){ 
	  setHomeDir(homeDir);
	  setDataDir(dataDir);
	  setTmpDir(tmpDir);
  }  
  
  public List<SeqRocket> launchSingleEnd(String inSeqFile,String libSeqInfoFile, 
			 SeqRocketConfig rocketConfig, String seqOutDir){	    
			
	  FileOperation.newFolder(tmpDir+"/recognizer");		
	  List<String> tmpFiles=new  ArrayList<String>();
					
	  if(seqOutDir==null){
		 seqOutDir=inSeqFile.substring(0, inSeqFile.lastIndexOf("/"))+"/RecognizedSeq";	
   	  }
      FileOperation.newFolder(seqOutDir);
	  seqOutDir=seqOutDir+"/";

	  // forward seq(Single-end)
	  List<SeqRocket> rockets =buildRockets(libSeqInfoFile,rocketConfig);		
					
	  List<SeqCompoAlignInfo> seqObjList=null; 
	  //List<SeqCompoAlignInfo> payLoadSeqObj;
					
	  //Setting and Checkpoint for pair-end forward Seq................
	  System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
	  seqObjList=checkSeq(inSeqFile);
	  System.out.println("Checked Seq Num: "+seqObjList.size());
	  String inSeqFileFormat=FileOperation.getFileFormat(inSeqFile);
			
	  //Recognizing barcode ...............................................		        
	  List<SeqCompoAlignInfo> barNonExactSeqObj
			       =getNonExactAndSaveExactSeq(seqObjList,rockets,BARCODE_NAME_DEFINITION,inSeqFileFormat);
	  tmpFiles.addAll(getCompoExactAlignedFiles(rockets,BARCODE_NAME_DEFINITION));
	  seqObjList=null;		
	  String barNonExactLeftSubSeqFile=tmpDir+"/AllBarNoExactMatchSeq_forward.leftsub.fna";
	  tmpFiles.add(barNonExactLeftSubSeqFile);
	  barcodeMaxTerritoryLen=getSeqCompoMaxTerritoryLen(rockets, BARCODE_NAME_DEFINITION);
	  createLeftSubBLASTTarSeq(barNonExactSeqObj,barcodeMaxTerritoryLen,barNonExactLeftSubSeqFile);	
			
	  for(int i=0;i<rockets.size();i++){			  
		 launchSeqRocket(rockets.get(i),barNonExactSeqObj,barNonExactLeftSubSeqFile,seqOutDir);
	  }// for Rocket
  	  barNonExactSeqObj=null;    
		    
	  //System.out.println("Delete temporary files");
      for(String tmpFile: tmpFiles){
	     FileOperation.delFile(tmpFile);
	  }
	  tmpFiles=null;    
      setSeqRockets(rockets);
		    
	  return rockets;
		   
  }


  public List<SeqRocket> splitLaunchSingleEnd(String inSeqFile,int splitStep,String libSeqInfoFile,
		 SeqRocketConfig rocketConfig,String splitedSeqOut, String combinedSeqOut){
	 
	 List<String> splitedSeqFiles=new ArrayList<String>();
	//Check seq format, and then split seq into multiple subfiles................
	 System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
	 splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile, splitStep, splitedSeqOut);	 
	 
	 List<SeqRocket> seqRockets=splitLaunchSingleEnd(splitedSeqFiles,libSeqInfoFile,rocketConfig,combinedSeqOut);
	 
	 return seqRockets;
  } 
 
  public List<SeqRocket> splitLaunchSingleEnd(List<String> splitedSeqFiles, String libSeqInfoFile, 
		 SeqRocketConfig rocketConfig, String combinedSeqOut){
     
	 if(combinedSeqOut==null) combinedSeqOut=GEAT.getWorkingDir(); 
	 
	 List<SeqRocket> seqRockets=new ArrayList<SeqRocket>();
	 
	 List<ArrayList<SeqRocket>> splitedRockets=new ArrayList<ArrayList<SeqRocket>>();
	 ArrayList<SeqRocket> rockets;
	 List<String> tmpFiles=new ArrayList<String>();
	 String outDir = null;
	 String splitDir=null;
	 List<Boolean> saveRecognizedSeq=new ArrayList<Boolean>();
	 List<Boolean> isRecognized=new ArrayList<Boolean>();
	 List<Boolean> isEndTrimmed=new ArrayList<Boolean>();
	 List<Boolean> isEndMasked=new ArrayList<Boolean>();
	 int s=1;
	 for(String file:splitedSeqFiles){
	    try {
		   splitDir=file.substring(0,file.lastIndexOf("/"));
		   outDir=combinedSeqOut+"/splits";
		   FileOperation.newFolder(outDir);
		   outDir=outDir+"/"+file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."));
		   FileOperation.newFolder(outDir);
		   
		   System.out.println("......Recognizing sequence for split "+s+" ......");
		   s++;
		   
		   rockets=(ArrayList<SeqRocket>) launchSingleEnd(file,libSeqInfoFile,rocketConfig,outDir);		  
	       if(rockets!=null && rockets.size()>0){
	    	  splitedRockets.add(rockets);
	    	  
	    	  for(SeqRocket rocket:rockets) {
	    		  saveRecognizedSeq.add(false);
	    		  isRecognized.add(false);
	    		  isEndTrimmed.add(false);
	    		  isEndMasked.add(false);
	    		  
				  if(rocket.savedCompoAlignedSeqFiles!=null) {
				     for(String savedCompoFile:rocket.savedCompoAlignedSeqFiles) {
				       tmpFiles.add(savedCompoFile);
				     }
				  }
				  if(rocket.recognizedSeqFile!=null)
				     tmpFiles.add(rocket.recognizedSeqFile); 
				  if(rocket.recognizedSeqFileMasked!=null)
				     tmpFiles.add(rocket.recognizedSeqFileMasked);
				  if(rocket.recognizedSeqFileTrimmed!=null)
				     tmpFiles.add(rocket.recognizedSeqFileTrimmed); 
	    	  }	    	  
		      tmpFiles.add(file);
		      
		      for(int r=0;r<rockets.size();r++){
		    	SeqRocket rocket=rockets.get(r);
		    	if(rocket.saveRecognizedSeq) saveRecognizedSeq.set(r, true);
		    	if(rocket.isRecognized) isRecognized.set(r, true);
		    	if(rocket.isEndTrimmed) isEndTrimmed.set(r, true);
		    	if(rocket.isEndMasked) isEndMasked.set(r, true);
			    rocket=null;
		      }
	       }
		   rockets=null;
		}catch (Exception e) {
		// TODO Auto-generated catch block
		   e.printStackTrace();
		}
	 }
	   
	 if(splitedSeqFiles.size()>0 && splitedRockets.size()>0){
		 String file=splitedSeqFiles.get(0);
		
		 FileOperation.newFolder(combinedSeqOut);
		 //System.out.println("......Combining splited sequences......");
		 List<String> fileList=new ArrayList<String>();
		 String outName="";
		 String seqFormat="fna";
	     for(int r=0;r<splitedRockets.get(0).size();r++){
	    	 SeqRocket rocket=new SeqRocket();
	    	 rocket=splitedRockets.get(0).get(r); 
			 rocket.saveRecognizedSeq=saveRecognizedSeq.get(r);
			 rocket.isRecognized=isRecognized.get(r);
			 rocket.isEndTrimmed=isEndTrimmed.get(r);
			 rocket.isEndMasked=isEndMasked.get(r);
		 
	    	 // for recognized Seq
	    	 if(rocket.isRecognized && rocket.saveRecognizedSeq){
		    	 System.out.println(rocket.rocketName+" >>> combining finally recognized seq......");
				 splitedSeqFiles=new ArrayList<String>();
				 for(int i=0;i<splitedRockets.size();i++){
					splitedSeqFiles.add(splitedRockets.get(i).get(r).recognizedSeqFile);
				 }
				 
				 file=rocket.rocketName+".";
				 for(String fileName:splitedSeqFiles) {
					if(fileName!=null && new File(fileName).exists()) {
					   file=fileName; 
					   break;
					}
				 }
				 seqFormat=FileOperation.getFileFormat(file);
				 outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles,outName);
				 rocket.recognizedSeqFile=outName;
				 fileList.add(outName);
	    	 }else {
	    		 if(rocket.isRecognized)
	   			   System.out.println(rocket.rocketName+" >>> Warning: you did not set to save finally recognized seq, please check '-finalSave'!!!");
	    		 else
		    	   System.out.println(rocket.rocketName+" >>> Warning: No finally recognized seq, or check if your data or configuration is possibly incorrect for it!!!");
	    	 }
			 
			 //for trimmed Seq file
			 if(rocket.isEndTrimmed){
			     System.out.println(rocket.rocketName+" >>> combining finally recognized seq trimmed......");
			     splitedSeqFiles=new ArrayList<String>();			
			     for(int i=0;i<splitedRockets.size();i++){
				    splitedSeqFiles.add(splitedRockets.get(i).get(r).recognizedSeqFileTrimmed);
			     }
			     file=rocket.rocketName+".Trimmed.";
				 for(String fileName:splitedSeqFiles) {
					if(fileName!=null && new File(fileName).exists()) {
					   file=fileName; 
					   break;
					}
				 }
			     seqFormat=FileOperation.getFileFormat(file);
			     outName=combinedSeqOut+"/"
			          +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		              +".final."+seqFormat;
			     SeqOperation.combineSeqFile(splitedSeqFiles,outName);
			     rocket.recognizedSeqFileTrimmed=outName;
			 }else {
	    		 if(rocket.isRecognized)
				   System.out.println(rocket.rocketName+" >>> Warning: you did not set to TRIM finally recognized file, please check the related parameters !!!");
			 }

			 
			 //for masked Seq file
			 if(rocket.isRecognized && rocket.isEndMasked){
			     System.out.println(rocket.rocketName+" >>> combining finally recognized seq masked......");
			     splitedSeqFiles=new ArrayList<String>();			
			     for(int i=0;i<splitedRockets.size();i++){
				   splitedSeqFiles.add(splitedRockets.get(i).get(r).recognizedSeqFileMasked);
			     }
			     file=rocket.rocketName+".Masked.";
				 for(String fileName:splitedSeqFiles) {
					if(fileName!=null && new File(fileName).exists()) {
					   file=fileName; 
					   break;
					}
				 }
			     seqFormat=FileOperation.getFileFormat(file);
			     outName=combinedSeqOut+"/"
			          +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		              +".final."+seqFormat;
			     SeqOperation.combineSeqFile(splitedSeqFiles,outName);
			     rocket.recognizedSeqFileMasked=outName;
			 }else {
	    		 if(rocket.isRecognized)
				   System.out.println(rocket.rocketName+" >>> Warning: you did not set to MASK finally recognized file, please check the related parameters !!!");
			 }
			 
			 //for saved component-recognized Seq
			 if(!rocket.finalSeqCompo.equalsIgnoreCase(BARCODE_NAME_DEFINITION)){
				for(int c=0;c<rocket.savedCompoAlignedSeqFiles.size();c++){
		    	  splitedSeqFiles=new ArrayList<String>();	    	 
				  for(int i=0;i<splitedRockets.size();i++){
					splitedSeqFiles.add(splitedRockets.get(i).get(r).savedCompoAlignedSeqFiles.get(c));
				  }
				  file=rocket.rocketName+".";
				  for(String fileName:splitedSeqFiles){
					if(fileName!=null && new File(fileName).exists()){
					   file=fileName; 
					   break;
					}
				  }
				  seqFormat=FileOperation.getFileFormat(file);
				  outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".required."+seqFormat;
				  SeqOperation.combineSeqFile(splitedSeqFiles,outName);
				  rocket.savedCompoAlignedSeqFiles.set(c,outName);
				}
			 }
			 
			 seqRockets.add(rocket);
			 rocket=null;
	     }
	     
	     String fileListOut=combinedSeqOut+"/fileList.txt";
		 FileOperation.saveList(fileList, null, fileListOut);

	 }else{
		System.out.println("Warning: No sequences recognized, please check if your data or configuration is correct!!!");
	 }
	 splitedRockets=null;
	 saveRecognizedSeq=null;
	 isRecognized=null;
	 isEndTrimmed=null;
	 isEndMasked=null;
	 
	 setSeqRockets(seqRockets);	
	 
	 if(tmpFiles!=null){
		for(String tmpFile: tmpFiles){
		   FileOperation.delFile(tmpFile);
		}
	 }
	 tmpFiles=null;
	 if(splitDir!=null) FileOperation.delFolder(splitDir);
	 
	 System.out.println("......Seq Recognization Done......");
	 
     return seqRockets;
  }
  
  boolean isRocketPairedMode(List<SeqRocket> rockets,List<SeqRocket> rockets2){
	  boolean isOK=true;
	  if(rockets.size()!=rockets2.size()) return false;
	  
	  for(int i=0;i<rockets.size();i++) {
		 if(!rockets.get(i).expName.equalsIgnoreCase(rockets2.get(i).expName)) {
			isOK=false;
			break;
		 }
	  }
	  return isOK;
  } 
  
  public List<SeqRocketPair> launchPairEnd(String inSeqFile,String inSeqFile2,
		   String libSeqInfoFile,String libSeqInfoFile2, String pairMode,
		   SeqRocketConfig rocketConfig, SeqRocketConfig rocketConfig2, 
		   String seqOutDir) {	
	  
	  List<SeqRocketPair> rockets12 = null;
	  
	  if(pairMode.equalsIgnoreCase(LIBSEQINFO_PEMODE_PW)){
	      rockets12=(ArrayList<SeqRocketPair>) launchPairEnd_PairWise(
	    	  inSeqFile,inSeqFile2,libSeqInfoFile,libSeqInfoFile2,rocketConfig,rocketConfig2,seqOutDir
			);
	  }else if(pairMode.equalsIgnoreCase(LIBSEQINFO_PEMODE_O2O)){
	      rockets12=(ArrayList<SeqRocketPair>) launchPairEnd_Pair(
	    	  inSeqFile,inSeqFile2,libSeqInfoFile,libSeqInfoFile2,rocketConfig,rocketConfig2,seqOutDir
			);
	  }
	  
	  return rockets12;
  }  
  
  List<SeqRocketPair> launchPairEnd_Pair(String inSeqFile,String inSeqFile2,
		   String libSeqInfoFile,String libSeqInfoFile2, 
		   SeqRocketConfig rocketConfig, SeqRocketConfig rocketConfig2, 
		   String seqOutDir) {	    
		
		FileOperation.newFolder(tmpDir+"/recognizer");		
		List<String>tmpFiles=new  ArrayList<String>();
				
		if(seqOutDir==null){
		   seqOutDir=inSeqFile.substring(0, inSeqFile.lastIndexOf("/"))+"/RecognizedSeq";	
		}
		FileOperation.newFolder(seqOutDir);
		seqOutDir=seqOutDir+"/";
				
		List<SeqCompoAlignInfo> seqObjList=null; //pair-end forward
		List<SeqCompoAlignInfo> seqObjList2=null; //pair-end reversed
		//List<SeqCompoAlignInfo> payLoadSeqObj;
				
		//Setting and Checkpoint for pair-end forward Seq................
		System.out.println("Total Seq Num(pair-end forward): "+SeqOperation.getSeqNum(inSeqFile));
		String inSeqFileFormat=FileOperation.getFileFormat(inSeqFile);
		seqObjList=checkSeq(inSeqFile);
		System.out.println("Checked Seq Num(pair-end forward): "+seqObjList.size());
		
		//Setting and Checkpoint for pair-end reversed Seq................
		System.out.println("Total Seq Num(pair-end reverse): "+SeqOperation.getSeqNum(inSeqFile2));
		String inSeqFileFormat2=FileOperation.getFileFormat(inSeqFile2);
		seqObjList2=checkSeq(inSeqFile2);
		System.out.println("Checked Seq Num(pair-end reverse): "+seqObjList2.size());
   
		System.out.println("Checking pair-end seq......");
		if(!checkPairEndSeq(seqObjList,seqObjList2)){
			System.err.println("The pair-end seq are not consistent each other.");
			return null;
		}		
		seqObjList2=null; // release memory space
		
		List<SeqRocketPair> rockets12=new ArrayList<SeqRocketPair>();// forward-reverse combination		
		List<SeqRocket> rockets =buildRockets(libSeqInfoFile,rocketConfig);	// for forward seq(Pair-end)
		List<SeqRocket> rockets2=buildRockets(libSeqInfoFile2,rocketConfig2); // for reversed seq(Pair-end)
		if(!isRocketPairedMode(rockets,rockets2)){
		   System.err.println("The pair-end experiment seq info of library are not consistent each other.");
		   return null;
		}
		
		//Recognizing barcode for pair-end forward Seq...............................................	
		List<SeqCompoAlignInfo> barNonExactSeqObj
		        =getNonExactAndSaveExactSeq(seqObjList,rockets,BARCODE_NAME_DEFINITION,inSeqFileFormat);
		tmpFiles.addAll(getCompoExactAlignedFiles(rockets,BARCODE_NAME_DEFINITION));
		seqObjList=null;
        String barNonExactLeftSubSeqFile=tmpDir+"/AllBarNoExactMatchSeq_forward.leftsub.fna";
        tmpFiles.add(barNonExactLeftSubSeqFile);
        barcodeMaxTerritoryLen=getSeqCompoMaxTerritoryLen(rockets, BARCODE_NAME_DEFINITION);
        createLeftSubBLASTTarSeq(barNonExactSeqObj,barcodeMaxTerritoryLen,barNonExactLeftSubSeqFile);       

		for(int i=0;i<rockets.size();i++){		  
		    
			SeqRocket rocket_forward=rockets.get(i);		    
            String rocket12Tag=rockets.get(i).expName;
		    rocket_forward.rocketName=rocket12Tag+"_PairEnd.Forward";
	  
			launchSeqRocket(rocket_forward,barNonExactSeqObj,barNonExactLeftSubSeqFile,seqOutDir);
			
			String forwardSeqFile=rockets.get(i).recognizedSeqFile;			
			if(forwardSeqFile==null || !new File(forwardSeqFile).exists()) {
			  // System.out.println("Warning: no finally saved or recognized forward file for "
			  //     +rockets.get(i).rocketName
			  //     +"!, please check if '-finalSave' or other input information is correct for it! We just skipped it!");
			   continue;
			}
			
			String forwardSeqFormat=FileOperation.getFileFormat(forwardSeqFile);
			String forwardSeqNameFile=tmpDir+"/"+forwardSeqFile.substring(
					 forwardSeqFile.lastIndexOf("/")+1,forwardSeqFile.lastIndexOf(".")
			    	)+".seqName";
		    tmpFiles.add(forwardSeqNameFile);
		    tmpFiles.add(forwardSeqFile);
		   
			for(String forwardSavedSeqFile:rockets.get(i).savedCompoAlignedSeqFiles) {		
			  if(!forwardSavedSeqFile.equals(forwardSeqFile)){		
			     tmpFiles.add(forwardSavedSeqFile);
			  }
			}
			 
			String forwardMaskSeqFile=rockets.get(i).recognizedSeqFileMasked;
			String forwardMaskSeqFormat=FileOperation.getFileFormat(forwardMaskSeqFile);
		    tmpFiles.add(forwardMaskSeqFile);
		   
			String forwardTrimSeqFile=rockets.get(i).recognizedSeqFileTrimmed;
			String forwardTrimSeqFormat=FileOperation.getFileFormat(forwardTrimSeqFile);
			tmpFiles.add(forwardTrimSeqFile);
			
			SeqOperation.extratSeqName(forwardSeqFile,forwardSeqNameFile);
			String reverseSeqFile=forwardSeqFile.substring(
					                0,forwardSeqFile.lastIndexOf(".")
					              )+".Reverse."+inSeqFileFormat2;
			SeqOperation.getSubSeq(inSeqFile2,forwardSeqNameFile,0,reverseSeqFile);
			seqObjList2=checkSeq(reverseSeqFile);
			tmpFiles.add(reverseSeqFile);			
					
	        //List<SeqRocket> rockets2=buildRockets(libSeqInfoFile2, rocketConfig2); // for reversed seq(Pair-end)
			List<SeqCompoAlignInfo>barNonExactSeqObj2
			    =getNonExactAndSaveExactSeq(seqObjList2,rockets2,BARCODE_NAME_DEFINITION,inSeqFileFormat2);
			tmpFiles.addAll(getCompoExactAlignedFiles(rockets2,BARCODE_NAME_DEFINITION));	
			seqObjList2=null;
			String barNonExactLeftSubSeqFile2=tmpDir+"/AllBarNoExactMatchSeq_reverse.leftsub.fna";
		    tmpFiles.add(barNonExactLeftSubSeqFile2);	
		    barcodeMaxTerritoryLen=getSeqCompoMaxTerritoryLen(rockets2, BARCODE_NAME_DEFINITION);
		    createLeftSubBLASTTarSeq(barNonExactSeqObj2,barcodeMaxTerritoryLen,barNonExactLeftSubSeqFile2);	
				    
			/////////////// get reverse-recognized result
			SeqRocket rocket_reverse=rockets2.get(i);
			rocket_reverse.rocketName=rocket12Tag+"_PairEnd";					
	        launchSeqRocket(rocket_reverse,barNonExactSeqObj2,barNonExactLeftSubSeqFile2,seqOutDir);
							
			////// get corresponding forward seq based on reverse-recognized result //////
			// for corresponding forward seq of finally reverse-recognized seq
		    String reverseSeqOut=rocket_reverse.recognizedSeqFile;	
		    String reverseSeqFormat=FileOperation.getFileFormat(reverseSeqOut);
			String reverseSeqNameFile = null;
			String reverseSeqOut2=null;		
			String forwardSeqOut = null;			
			if(reverseSeqOut!=null && new File(reverseSeqOut).exists()) {
			   reverseSeqNameFile=tmpDir+"/"+reverseSeqOut.substring(
						  reverseSeqOut.lastIndexOf("/")+1,reverseSeqOut.lastIndexOf(".")
				     )+".seqName";	
			   tmpFiles.add(reverseSeqNameFile);			       
			   SeqOperation.extratSeqName(reverseSeqOut,reverseSeqNameFile);					
			   forwardSeqOut=reverseSeqOut.substring(
					      0,reverseSeqOut.lastIndexOf(".")
					 )+".Forward."+forwardSeqFormat;				
			   SeqOperation.getSubSeq(forwardSeqFile,reverseSeqNameFile,0,forwardSeqOut);
			   rocket_forward.recognizedSeqFile=forwardSeqOut;
			   
			   reverseSeqOut2=reverseSeqOut.substring(
					      0,reverseSeqOut.lastIndexOf(".")
					 )+".Reverse."+reverseSeqFormat;
			   new File(reverseSeqOut).renameTo(new File(reverseSeqOut2));
			   rocket_reverse.recognizedSeqFile=reverseSeqOut2;
			   reverseSeqOut2=null;
			}else {
			   rocket_forward.recognizedSeqFile=null;
			   rocket_reverse.recognizedSeqFile=null;
			   if(rocket_forward.isRecognized) rocket_forward.saveRecognizedSeq=false;
			   if(rocket_reverse.isRecognized) rocket_reverse.saveRecognizedSeq=false;
			}
					
			// for corresponding forward seq of reverse component-recognized seq to be saved				  
			rocket_forward.savedCompoAlignedSeqFiles=new ArrayList<String>();
			for(int s=0;s<rocket_reverse.savedCompoAlignedSeqFiles.size();s++) {
			   String revSavedSeqOut=rocket_reverse.savedCompoAlignedSeqFiles.get(s);			
			   String revSavedSeqFormat=FileOperation.getFileFormat(revSavedSeqOut);
			   String revSavedSeqOut2=null;
			   
			   String forSavedSeqFile=forwardSeqFile;
			   String forSavedSeqFormat=forwardSeqFormat;
			   
			   if(revSavedSeqOut.equals(reverseSeqOut)){
				  rocket_forward.savedCompoAlignedSeqFiles.add(forwardSeqOut);
			   }else{	
				  String revSavedSeqNameFile=tmpDir+"/"+revSavedSeqOut.substring(
						  revSavedSeqOut.lastIndexOf("/")+1,revSavedSeqOut.lastIndexOf(".")
				      )+".seqName";	
				  tmpFiles.add(revSavedSeqNameFile);			       
				  SeqOperation.extratSeqName(revSavedSeqOut,revSavedSeqNameFile);					
				  String forSavedSeqOut=revSavedSeqOut.substring(
				          0,revSavedSeqOut.lastIndexOf(".")
				       )+".Forward."+forSavedSeqFormat;				
				  SeqOperation.getSubSeq(forSavedSeqFile,revSavedSeqNameFile,0,forSavedSeqOut);
				  rocket_forward.savedCompoAlignedSeqFiles.add(forSavedSeqOut);
				  forSavedSeqOut=null;
				  revSavedSeqNameFile=null;
			   }
			   forSavedSeqFile=null;
			   forSavedSeqFormat=null;
			   
			   revSavedSeqOut2=revSavedSeqOut.substring(
					      0,revSavedSeqOut.lastIndexOf(".")
					 )+".Reverse."+revSavedSeqFormat;
			   new File(revSavedSeqOut).renameTo(new File(revSavedSeqOut2));
			   rocket_reverse.savedCompoAlignedSeqFiles.set(s,revSavedSeqOut2); 
			   revSavedSeqOut2=null;
			   revSavedSeqOut=null;	
			   revSavedSeqFormat=null;
			}
					
			// for corresponding forward trimmed seq of finally reverse-recognized seq
			if(rocket_forward.isEndTrimmed && rocket_reverse.isEndTrimmed) {
				String reverseTrimSeqOut=rocket_reverse.recognizedSeqFileTrimmed;
			    String reverseTrimSeqFormat=FileOperation.getFileFormat(reverseTrimSeqOut);
				String reverseTrimSeqOut2;
				 
				if(reverseSeqNameFile==null || !new File(reverseSeqNameFile).exists()) {				  
				   reverseSeqNameFile=tmpDir+"/"+reverseTrimSeqOut.substring(
						   reverseTrimSeqOut.lastIndexOf("/")+1,reverseTrimSeqOut.lastIndexOf(".")
				      )+".seqName";	
			       if(new File(reverseTrimSeqOut).exists()) {
					  tmpFiles.add(reverseSeqNameFile);				      
					  SeqOperation.extratSeqName(reverseTrimSeqOut,reverseSeqNameFile);
				   }
				}
				
				String forwardTrimSeqOut=rocket_reverse.recognizedSeqFileTrimmed.substring(
					   0,rocket_reverse.recognizedSeqFileTrimmed.lastIndexOf(".")
				    )+".Forward."+forwardTrimSeqFormat;
			    SeqOperation.getSubSeq(forwardTrimSeqFile,reverseSeqNameFile,0,forwardTrimSeqOut);
			    rocket_forward.recognizedSeqFileTrimmed=forwardTrimSeqOut;	
			    forwardTrimSeqOut=null;
			     
			    reverseTrimSeqOut2=reverseTrimSeqOut.substring(
				       0,reverseTrimSeqOut.lastIndexOf(".")
					 )+".Reverse."+reverseTrimSeqFormat;
				new File(reverseTrimSeqOut).renameTo(new File(reverseTrimSeqOut2));
				rocket_reverse.recognizedSeqFileTrimmed=reverseTrimSeqOut2;
				reverseTrimSeqOut2=null;
				reverseTrimSeqFormat=null;
			}else {
			    rocket_forward.recognizedSeqFileTrimmed=null;
			    rocket_reverse.recognizedSeqFileTrimmed=null;
			}
					
		    // for corresponding forward masked seq of finally reverse-recognized seq
			if(rocket_forward.isEndMasked && rocket_reverse.isEndMasked) {
			    String reverseMaskSeqOut=rocket_reverse.recognizedSeqFileMasked;
				String reverseMaskSeqFormat=FileOperation.getFileFormat(reverseMaskSeqOut);
				String reverseMaskSeqOut2;
			    if(reverseSeqNameFile==null || !new File(reverseSeqNameFile).exists()) {				
					reverseSeqNameFile=tmpDir+"/"+reverseMaskSeqOut.substring(
							reverseMaskSeqOut.lastIndexOf("/")+1,reverseMaskSeqOut.lastIndexOf(".")
					      )+".seqName";	
					if(new File(reverseMaskSeqOut).exists()) {
					   tmpFiles.add(reverseSeqNameFile);				      
					   SeqOperation.extratSeqName(reverseMaskSeqOut,reverseSeqNameFile);
					}
			    }
				String forwardMaskSeqOut=rocket_reverse.recognizedSeqFileMasked.substring(
				       0,rocket_reverse.recognizedSeqFileMasked.lastIndexOf(".")
				    )+".Forward."+forwardMaskSeqFormat;
				SeqOperation.getSubSeq(forwardMaskSeqFile,reverseSeqNameFile,0,forwardMaskSeqOut);
				rocket_forward.recognizedSeqFileMasked=forwardMaskSeqOut;	
				forwardMaskSeqOut=null;
				
				reverseMaskSeqOut2=reverseMaskSeqOut.substring(
					   0,reverseMaskSeqOut.lastIndexOf(".")
					)+".Reverse."+reverseMaskSeqFormat;
			    new File(reverseMaskSeqOut).renameTo(new File(reverseMaskSeqOut2));
			    rocket_reverse.recognizedSeqFileMasked=reverseMaskSeqOut2;
			    reverseMaskSeqOut2=null;
			    reverseMaskSeqFormat=null;
			}else {
			    rocket_forward.recognizedSeqFileMasked=null;
			    rocket_reverse.recognizedSeqFileMasked=null;
			}
					
			SeqRocketPair rocketPair=new SeqRocketPair();
			rocket_forward.rocketName=rocket12Tag+"_PairEnd.Forward";
			rocket_reverse.rocketName=rocket12Tag+"_PairEnd.Reverse";
			rocketPair.name=rocket12Tag;
			rocketPair.forward=rocket_forward;
			rocketPair.reverse=rocket_reverse;
			rockets12.add(rocketPair);
					
			rocketPair=null;
			rocket12Tag=null;
			reverseSeqOut=null;
			reverseSeqFormat=null;
			reverseSeqNameFile=null;									
			forwardSeqOut=null;	
			
			barNonExactSeqObj2=null;
			barNonExactLeftSubSeqFile2=null;
			reverseSeqFile=null;
			forwardSeqFile=null;
			forwardSeqFormat=null;
			forwardSeqNameFile=null;			
			forwardMaskSeqFile=null;
			forwardMaskSeqFormat=null;
			forwardTrimSeqFile=null;
			forwardTrimSeqFormat=null;		
					
		}// for Rocket
		barNonExactSeqObj=null;    
		barNonExactLeftSubSeqFile=null;
		rockets=null;
		rockets2=null;
	    setSeqPairRockets(rockets12);	    
	 
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
	      FileOperation.delFile(tmpFile);
		}
	    tmpFiles=null;
	    
	    return rockets12;
	  
  }
  
  List<SeqRocketPair> launchPairEnd_PairWise(String inSeqFile,String inSeqFile2,
		   String libSeqInfoFile,String libSeqInfoFile2, 
		   SeqRocketConfig rocketConfig, SeqRocketConfig rocketConfig2, 
		   String seqOutDir) {	    
		
		FileOperation.newFolder(tmpDir+"/recognizer");		
		List<String>tmpFiles=new  ArrayList<String>();
				
		if(seqOutDir==null){
		   seqOutDir=inSeqFile.substring(0, inSeqFile.lastIndexOf("/"))+"/RecognizedSeq";	
		}
		FileOperation.newFolder(seqOutDir);
		seqOutDir=seqOutDir+"/";
				
		List<SeqCompoAlignInfo> seqObjList=null; //pair-end forward
		List<SeqCompoAlignInfo> seqObjList2=null; //pair-end reversed
		//List<SeqCompoAlignInfo> payLoadSeqObj;
				
		//Setting and Checkpoint for pair-end forward Seq................
		System.out.println("Total Seq Num(pair-end forward): "+SeqOperation.getSeqNum(inSeqFile));
		String inSeqFileFormat=FileOperation.getFileFormat(inSeqFile);
		seqObjList=checkSeq(inSeqFile);
		System.out.println("Checked Seq Num(pair-end forward): "+seqObjList.size());
		
		// Setting and Checkpoint for pair-end reversed Seq................
		System.out.println("Total Seq Num(pair-end reverse): "+SeqOperation.getSeqNum(inSeqFile2));
		String inSeqFileFormat2=FileOperation.getFileFormat(inSeqFile2);
		seqObjList2=checkSeq(inSeqFile2);
		System.out.println("Checked Seq Num(pair-end reverse): "+seqObjList2.size());
      ///*
		System.out.println("Checking pair-end seq......");
		if(!checkPairEndSeq(seqObjList,seqObjList2)){
			System.err.println("The pair-end seq are not consistent each other.");
			return null;
		}
		//*/
		seqObjList2=null; // release memory space
		
		List<SeqRocketPair> rockets12=new ArrayList<SeqRocketPair>();// forward-reverse combination		
		List<SeqRocket> rockets =buildRockets(libSeqInfoFile,rocketConfig);	// for forward seq(Pair-end)
		
		//Recognizing barcode for pair-end forward Seq...............................................	
		List<SeqCompoAlignInfo> barNonExactSeqObj
		        =getNonExactAndSaveExactSeq(seqObjList,rockets,BARCODE_NAME_DEFINITION,inSeqFileFormat);
		tmpFiles.addAll(getCompoExactAlignedFiles(rockets,BARCODE_NAME_DEFINITION));
		seqObjList=null;
       String barNonExactLeftSubSeqFile=tmpDir+"/AllBarNoExactMatchSeq_forward.leftsub.fna";
       tmpFiles.add(barNonExactLeftSubSeqFile);
       barcodeMaxTerritoryLen=getSeqCompoMaxTerritoryLen(rockets, BARCODE_NAME_DEFINITION);
       createLeftSubBLASTTarSeq(barNonExactSeqObj,barcodeMaxTerritoryLen,barNonExactLeftSubSeqFile);	

		for(int i=0;i<rockets.size();i++){		  
		  
			launchSeqRocket(rockets.get(i),barNonExactSeqObj,barNonExactLeftSubSeqFile,seqOutDir);
			
			String forwardSeqFile=rockets.get(i).recognizedSeqFile;
			
			if(forwardSeqFile==null || !new File(forwardSeqFile).exists()) {
			  // System.out.println("Warning: no finally saved or recognized forward file for "
			  //     +rockets.get(i).rocketName
			  //     +"!, please check if '-finalSave' or other input information is correct for it! We just skipped it!");
			   continue;
			}
			
			String forwardSeqFormat=FileOperation.getFileFormat(forwardSeqFile);
			String forwardSeqNameFile=tmpDir+"/"+forwardSeqFile.substring(
					 forwardSeqFile.lastIndexOf("/")+1,forwardSeqFile.lastIndexOf(".")
			    	)+".seqName";
		    tmpFiles.add(forwardSeqNameFile);
		    tmpFiles.add(forwardSeqFile);
		   
			for(String forwardSavedSeqFile:rockets.get(i).savedCompoAlignedSeqFiles) {		
			  if(!forwardSavedSeqFile.equals(forwardSeqFile)){		
			     tmpFiles.add(forwardSavedSeqFile);
			  }
			}
			 
			String forwardMaskSeqFile=rockets.get(i).recognizedSeqFileMasked;
			String forwardMaskSeqFormat=FileOperation.getFileFormat(forwardMaskSeqFile);
		    tmpFiles.add(forwardMaskSeqFile);
		   
			String forwardTrimSeqFile=rockets.get(i).recognizedSeqFileTrimmed;
			String forwardTrimSeqFormat=FileOperation.getFileFormat(forwardTrimSeqFile);
			tmpFiles.add(forwardTrimSeqFile);
			
			SeqOperation.extratSeqName(forwardSeqFile,forwardSeqNameFile);
			String reverseSeqFile=forwardSeqFile.substring(
					                0,forwardSeqFile.lastIndexOf(".")
					              )+".Reverse."+inSeqFileFormat2;
			SeqOperation.getSubSeq(inSeqFile2,forwardSeqNameFile,0,reverseSeqFile);
			seqObjList2=checkSeq(reverseSeqFile);
			tmpFiles.add(reverseSeqFile);			
					
	        List<SeqRocket> rockets2=buildRockets(libSeqInfoFile2, rocketConfig2); // for reversed seq(Pair-end)
			List<SeqCompoAlignInfo>barNonExactSeqObj2
			    =getNonExactAndSaveExactSeq(seqObjList2,rockets2,BARCODE_NAME_DEFINITION,inSeqFileFormat2);
			tmpFiles.addAll(getCompoExactAlignedFiles(rockets2,BARCODE_NAME_DEFINITION));	
			seqObjList2=null;
			String barNonExactLeftSubSeqFile2=tmpDir+"/AllBarNoExactMatchSeq_reverse.leftsub.fna";
		    tmpFiles.add(barNonExactLeftSubSeqFile2);	
		    barcodeMaxTerritoryLen=getSeqCompoMaxTerritoryLen(rockets2, BARCODE_NAME_DEFINITION);
		    createLeftSubBLASTTarSeq(barNonExactSeqObj2,barcodeMaxTerritoryLen,barNonExactLeftSubSeqFile2);	

			for(int j=0;j<rockets2.size();j++){	
			   try {  
				    
				    String rocket12Tag=rockets.get(i).rocketName+"-"+rockets2.get(j).rocketName;
				  
				    // get corresponding forward seq based on reverse-recognized result //////
				    SeqRocket rocket_forward=(SeqRocket) rockets.get(i).clone();
				  
				    // get reverse-recognized result
				    SeqRocket rocket_reverse=rockets2.get(j);
				    rocket_reverse.rocketName=rocket12Tag+"_PairEnd";					
		            launchSeqRocket(rocket_reverse,barNonExactSeqObj2,barNonExactLeftSubSeqFile2,seqOutDir);
					
					////// get corresponding forward seq based on reverse-recognized result //////
					// for corresponding forward seq of finally reverse-recognized seq
				    String reverseSeqOut=rocket_reverse.recognizedSeqFile;	
				    String reverseSeqFormat=FileOperation.getFileFormat(reverseSeqOut);
					String reverseSeqNameFile = null;
					String reverseSeqOut2=null;		
					String forwardSeqOut = null;			
					if(reverseSeqOut!=null && new File(reverseSeqOut).exists()) {
					   reverseSeqNameFile=tmpDir+"/"+reverseSeqOut.substring(
								  reverseSeqOut.lastIndexOf("/")+1,reverseSeqOut.lastIndexOf(".")
						     )+".seqName";	
					   tmpFiles.add(reverseSeqNameFile);			       
					   SeqOperation.extratSeqName(reverseSeqOut,reverseSeqNameFile);					
					   forwardSeqOut=reverseSeqOut.substring(
							      0,reverseSeqOut.lastIndexOf(".")
							 )+".Forward."+forwardSeqFormat;				
					   SeqOperation.getSubSeq(forwardSeqFile,reverseSeqNameFile,0,forwardSeqOut);
					   rocket_forward.recognizedSeqFile=forwardSeqOut;
					   
					   reverseSeqOut2=reverseSeqOut.substring(
							      0,reverseSeqOut.lastIndexOf(".")
							 )+".Reverse."+reverseSeqFormat;
					   new File(reverseSeqOut).renameTo(new File(reverseSeqOut2));
					   rocket_reverse.recognizedSeqFile=reverseSeqOut2;
					   reverseSeqOut2=null;
					}else {
					   rocket_forward.recognizedSeqFile=null;
					   rocket_reverse.recognizedSeqFile=null;
					   if(rocket_forward.isRecognized) rocket_forward.saveRecognizedSeq=false;
					   if(rocket_reverse.isRecognized) rocket_reverse.saveRecognizedSeq=false;
					}
							
					// for corresponding forward seq of reverse component-recognized seq to be saved				  
					rocket_forward.savedCompoAlignedSeqFiles=new ArrayList<String>();
					for(int s=0;s<rocket_reverse.savedCompoAlignedSeqFiles.size();s++) {
					   String revSavedSeqOut=rocket_reverse.savedCompoAlignedSeqFiles.get(s);			
					   String revSavedSeqFormat=FileOperation.getFileFormat(revSavedSeqOut);
					   String revSavedSeqOut2=null;
					   
					   String forSavedSeqFile=forwardSeqFile;
					   String forSavedSeqFormat=forwardSeqFormat;
					   
					   if(revSavedSeqOut.equals(reverseSeqOut)){
						  rocket_forward.savedCompoAlignedSeqFiles.add(forwardSeqOut);
					   }else{	
						  String revSavedSeqNameFile=tmpDir+"/"+revSavedSeqOut.substring(
								  revSavedSeqOut.lastIndexOf("/")+1,revSavedSeqOut.lastIndexOf(".")
						      )+".seqName";	
						  tmpFiles.add(revSavedSeqNameFile);			       
						  SeqOperation.extratSeqName(revSavedSeqOut,revSavedSeqNameFile);					
						  String forSavedSeqOut=revSavedSeqOut.substring(
						          0,revSavedSeqOut.lastIndexOf(".")
						       )+".Forward."+forSavedSeqFormat;				
						  SeqOperation.getSubSeq(forSavedSeqFile,revSavedSeqNameFile,0,forSavedSeqOut);
						  rocket_forward.savedCompoAlignedSeqFiles.add(forSavedSeqOut);
						  forSavedSeqOut=null;
						  revSavedSeqNameFile=null;
					   }
					   forSavedSeqFile=null;
					   forSavedSeqFormat=null;
					   
					   revSavedSeqOut2=revSavedSeqOut.substring(
							      0,revSavedSeqOut.lastIndexOf(".")
							 )+".Reverse."+revSavedSeqFormat;
					   new File(revSavedSeqOut).renameTo(new File(revSavedSeqOut2));
					   rocket_reverse.savedCompoAlignedSeqFiles.set(s,revSavedSeqOut2); 
					   revSavedSeqOut2=null;
					   revSavedSeqOut=null;	
					   revSavedSeqFormat=null;
					}
							
					// for corresponding forward trimmed seq of finally reverse-recognized seq
					if(rocket_forward.isEndTrimmed && rocket_reverse.isEndTrimmed) {
						String reverseTrimSeqOut=rocket_reverse.recognizedSeqFileTrimmed;
					    String reverseTrimSeqFormat=FileOperation.getFileFormat(reverseTrimSeqOut);
						String reverseTrimSeqOut2;
						 
						if(reverseSeqNameFile==null || !new File(reverseSeqNameFile).exists()) {				  
						   reverseSeqNameFile=tmpDir+"/"+reverseTrimSeqOut.substring(
								   reverseTrimSeqOut.lastIndexOf("/")+1,reverseTrimSeqOut.lastIndexOf(".")
						      )+".seqName";	
					       if(new File(reverseTrimSeqOut).exists()) {
							  tmpFiles.add(reverseSeqNameFile);				      
							  SeqOperation.extratSeqName(reverseTrimSeqOut,reverseSeqNameFile);
						   }
						}
						
						String forwardTrimSeqOut=rocket_reverse.recognizedSeqFileTrimmed.substring(
							   0,rocket_reverse.recognizedSeqFileTrimmed.lastIndexOf(".")
						    )+".Forward."+forwardTrimSeqFormat;
					    SeqOperation.getSubSeq(forwardTrimSeqFile,reverseSeqNameFile,0,forwardTrimSeqOut);
					    rocket_forward.recognizedSeqFileTrimmed=forwardTrimSeqOut;	
					    forwardTrimSeqOut=null;
					     
					    reverseTrimSeqOut2=reverseTrimSeqOut.substring(
						       0,reverseTrimSeqOut.lastIndexOf(".")
							 )+".Reverse."+reverseTrimSeqFormat;
						new File(reverseTrimSeqOut).renameTo(new File(reverseTrimSeqOut2));
						rocket_reverse.recognizedSeqFileTrimmed=reverseTrimSeqOut2;
						reverseTrimSeqOut2=null;
						reverseTrimSeqFormat=null;
					}else {
					    rocket_forward.recognizedSeqFileTrimmed=null;
					    rocket_reverse.recognizedSeqFileTrimmed=null;
					}
							
				    // for corresponding forward masked seq of finally reverse-recognized seq
					if(rocket_forward.isEndMasked && rocket_reverse.isEndMasked) {
					    String reverseMaskSeqOut=rocket_reverse.recognizedSeqFileMasked;
						String reverseMaskSeqFormat=FileOperation.getFileFormat(reverseMaskSeqOut);
						String reverseMaskSeqOut2;
					    if(reverseSeqNameFile==null || !new File(reverseSeqNameFile).exists()) {				
							reverseSeqNameFile=tmpDir+"/"+reverseMaskSeqOut.substring(
									reverseMaskSeqOut.lastIndexOf("/")+1,reverseMaskSeqOut.lastIndexOf(".")
							      )+".seqName";	
							if(new File(reverseMaskSeqOut).exists()) {
							   tmpFiles.add(reverseSeqNameFile);				      
							   SeqOperation.extratSeqName(reverseMaskSeqOut,reverseSeqNameFile);
							}
					    }
						String forwardMaskSeqOut=rocket_reverse.recognizedSeqFileMasked.substring(
						       0,rocket_reverse.recognizedSeqFileMasked.lastIndexOf(".")
						    )+".Forward."+forwardMaskSeqFormat;
						SeqOperation.getSubSeq(forwardMaskSeqFile,reverseSeqNameFile,0,forwardMaskSeqOut);
						rocket_forward.recognizedSeqFileMasked=forwardMaskSeqOut;	
						forwardMaskSeqOut=null;
						
						reverseMaskSeqOut2=reverseMaskSeqOut.substring(
							   0,reverseMaskSeqOut.lastIndexOf(".")
							)+".Reverse."+reverseMaskSeqFormat;
					    new File(reverseMaskSeqOut).renameTo(new File(reverseMaskSeqOut2));
					    rocket_reverse.recognizedSeqFileMasked=reverseMaskSeqOut2;
					    reverseMaskSeqOut2=null;
					    reverseMaskSeqFormat=null;
					}else {
					    rocket_forward.recognizedSeqFileMasked=null;
					    rocket_reverse.recognizedSeqFileMasked=null;
					}
							
					SeqRocketPair rocketPair=new SeqRocketPair();
					rocket_forward.rocketName=rocket12Tag+"_PairEnd.Forward";
					rocket_reverse.rocketName=rocket12Tag+"_PairEnd.Reverse";
					rocketPair.name=rocket12Tag;
					rocketPair.forward=rocket_forward;
					rocketPair.reverse=rocket_reverse;
					rockets12.add(rocketPair);
					
				    rocketPair=null;
				    rocket12Tag=null;
				    reverseSeqOut=null;
				    reverseSeqFormat=null;
				    reverseSeqNameFile=null;									
				    forwardSeqOut=null; 
					
			   }catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
				  e.printStackTrace();
			   }			
			}
			rockets2=null;
			barNonExactSeqObj2=null;
			barNonExactLeftSubSeqFile2=null;
			reverseSeqFile=null;
			forwardSeqFile=null;
			forwardSeqFormat=null;
			forwardSeqNameFile=null;			
			forwardMaskSeqFile=null;
			forwardMaskSeqFormat=null;
			forwardTrimSeqFile=null;
			forwardTrimSeqFormat=null;		
					
		}// for Rocket
		barNonExactSeqObj=null;    
		barNonExactLeftSubSeqFile=null;

	    setSeqPairRockets(rockets12);	    
	 
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
	      FileOperation.delFile(tmpFile);
		}
	    tmpFiles=null;
	    
	    return rockets12;
	  
  }
 
  public List<SeqRocketPair> splitLaunchPairEnd(String inSeqFile,String inSeqFile2,
		 String libSeqInfoFile,String libSeqInfoFile2, String pairMode,
		 SeqRocketConfig rocketConfig, SeqRocketConfig rocketConfig2,
		 int splitStep, String splitedSeqOut,String combinedSeqOut){
	 
	 List<String> splitedSeqFiles=new ArrayList<String>();
	 List<String> splitedSeqFiles2=new ArrayList<String>();
	 
	//Check seq format, and then split seq into multiple subfiles................	 
	 System.out.println("Total Seq Num(Pair-end Forward): "+SeqOperation.getSeqNum(inSeqFile));
	 splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile, splitStep, splitedSeqOut);
	 System.out.println("Total Seq Num(Pair-end Reverse): "+SeqOperation.getSeqNum(inSeqFile2));
	 splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2, splitStep, splitedSeqOut);

	 List<SeqRocketPair>seqPairRockets=splitLaunchPairEnd(splitedSeqFiles,splitedSeqFiles2,
			 libSeqInfoFile,libSeqInfoFile2,pairMode,rocketConfig, rocketConfig2,combinedSeqOut);	
	 
	 return seqPairRockets;

  }
 
  public List<SeqRocketPair> splitLaunchPairEnd(List<String>splitedSeqFiles, List<String>splitedSeqFiles2,
		 String libSeqInfoFile, String libSeqInfoFile2, String pairMode,
		 SeqRocketConfig rocketConfig, SeqRocketConfig rocketConfig2,
		 String combinedSeqOut){	 
     
	 if(combinedSeqOut==null) combinedSeqOut=GEAT.getWorkingDir();
	 List<SeqRocketPair> seqPairRockets=new ArrayList<SeqRocketPair>();
	 List<ArrayList<SeqRocketPair>> splitedRockets12=new ArrayList<ArrayList<SeqRocketPair>>();		 
	 ArrayList<SeqRocketPair> rockets12 = null;
	 List<String> tmpFiles=new ArrayList<String>();
	 String outDir =null;	
	 String file = null;
	 String file2 = null;
	 String forSplitDir = null;
	 String revSplitDir = null;
	 List<Boolean> saveRecognizedSeq=new ArrayList<Boolean>();
	 List<Boolean> isRecognized=new ArrayList<Boolean>();
	 List<Boolean> isEndTrimmed=new ArrayList<Boolean>();
	 List<Boolean> isEndMasked=new ArrayList<Boolean>();
	 List<Boolean> saveRecognizedSeq2=new ArrayList<Boolean>();
	 List<Boolean> isRecognized2=new ArrayList<Boolean>();
	 List<Boolean> isEndTrimmed2=new ArrayList<Boolean>();
	 List<Boolean> isEndMasked2=new ArrayList<Boolean>();
	 
	 for(int s=0;s<splitedSeqFiles.size();s++){
	    try {
	       file=splitedSeqFiles.get(s);
	       file2=splitedSeqFiles2.get(s);
	       
	       forSplitDir=file.substring(0,file.lastIndexOf("/"));
	       revSplitDir=file2.substring(0,file2.lastIndexOf("/"));
	       
	       outDir = combinedSeqOut+"/splits";
	       FileOperation.newFolder(outDir);
		   outDir=outDir+"/"+file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."));
		   FileOperation.newFolder(outDir);
		   
		   System.out.println("......Recognizing sequence for split "+(s+1)+" ......");
		   rockets12=(ArrayList<SeqRocketPair>) launchPairEnd(
			     file,file2,libSeqInfoFile,libSeqInfoFile2,pairMode,rocketConfig,rocketConfig2,outDir
			   );

		   if(rockets12!=null && rockets12.size()>0){
			   splitedRockets12.add(rockets12);		   
			   tmpFiles.add(file);
			   tmpFiles.add(file2);
			   for(SeqRocketPair rocketPair:rockets12){
				   
		    	  saveRecognizedSeq.add(false);
		    	  isRecognized.add(false);
		    	  isEndTrimmed.add(false);
		    	  isEndMasked.add(false);
		    		  
		    	  saveRecognizedSeq2.add(false);
		    	  isRecognized2.add(false);
		    	  isEndTrimmed2.add(false);
		    	  isEndMasked2.add(false);

                  // for forward
				  if(rocketPair.forward.savedCompoAlignedSeqFiles!=null) {
				     for(String savedCompoFile:rocketPair.forward.savedCompoAlignedSeqFiles) {
				       tmpFiles.add(savedCompoFile);
				     }
				  }
				  if(rocketPair.forward.recognizedSeqFile!=null)
				    tmpFiles.add(rocketPair.forward.recognizedSeqFile); 
				  if(rocketPair.forward.recognizedSeqFileMasked!=null)
				    tmpFiles.add(rocketPair.forward.recognizedSeqFileMasked);
				  if(rocketPair.forward.recognizedSeqFileTrimmed!=null)
				    tmpFiles.add(rocketPair.forward.recognizedSeqFileTrimmed); 
				  
                  // for reverse
				  if(rocketPair.reverse.savedCompoAlignedSeqFiles!=null) {
				    for(String savedCompoFile:rocketPair.reverse.savedCompoAlignedSeqFiles) {
				      tmpFiles.add(savedCompoFile);
				    }
				  }
				  if(rocketPair.reverse.recognizedSeqFile!=null)
					tmpFiles.add(rocketPair.reverse.recognizedSeqFile); 
				  if(rocketPair.reverse.recognizedSeqFileMasked!=null)
				    tmpFiles.add(rocketPair.reverse.recognizedSeqFileMasked);
				  if(rocketPair.reverse.recognizedSeqFileTrimmed!=null)
				    tmpFiles.add(rocketPair.reverse.recognizedSeqFileTrimmed); 
			   }
			   
			   for(int r=0;r<rockets12.size();r++){
				  SeqRocket rocket=rockets12.get(r).forward;
				  if(rocket.saveRecognizedSeq) saveRecognizedSeq.set(r,true);
				  if(rocket.isRecognized) isRecognized.set(r,true);
				  if(rocket.isEndTrimmed) isEndTrimmed.set(r,true);
				  if(rocket.isEndMasked) isEndMasked.set(r,true);
			      rocket=null;
			      SeqRocket rocket2=rockets12.get(r).reverse;
				  if(rocket2.saveRecognizedSeq) saveRecognizedSeq2.set(r,true);
				  if(rocket2.isRecognized) isRecognized2.set(r,true);
				  if(rocket2.isEndTrimmed) isEndTrimmed2.set(r,true);
				  if(rocket2.isEndMasked) isEndMasked2.set(r,true);
			      rocket2=null;
			   }

		   }
		   rockets12=null;
	
		}catch (Exception e) {
		// TODO Auto-generated catch block
		   e.printStackTrace();
		}
	 }
	   
	 if(splitedRockets12.size()>0){			
		 FileOperation.newFolder(combinedSeqOut);
		 System.out.println("......Combine splited sequences......");
		 List<String> fileList_forward=new ArrayList<String>();	
		 List<String> fileList_reverse=new ArrayList<String>();		
		 String outName="";
		 String outName2="";
		 String seqFormat="fna";
		 for(int r=0;r<splitedRockets12.get(0).size();r++){
			 SeqRocketPair rocketPair=new SeqRocketPair();
			 rocketPair.name=splitedRockets12.get(0).get(r).name;
			 rocketPair.forward=splitedRockets12.get(0).get(r).forward;			 
			 rocketPair.forward.saveRecognizedSeq=saveRecognizedSeq.get(r);
			 rocketPair.forward.isRecognized=isRecognized.get(r);
			 rocketPair.forward.isEndTrimmed=isEndTrimmed.get(r);
			 rocketPair.forward.isEndMasked=isEndMasked.get(r);
			 
			 rocketPair.reverse=splitedRockets12.get(0).get(r).reverse;
			 rocketPair.reverse.saveRecognizedSeq=saveRecognizedSeq2.get(r);
			 rocketPair.reverse.isRecognized=isRecognized2.get(r);
			 rocketPair.reverse.isEndTrimmed=isEndTrimmed2.get(r);
			 rocketPair.reverse.isEndMasked=isEndMasked2.get(r);			 

			 List<String> splitedSeqFiles_F;
			 List<String> splitedSeqFiles_R;
			 // for recognized Seq
			 if(rocketPair.forward.isRecognized && rocketPair.reverse.isRecognized
					&& rocketPair.forward.saveRecognizedSeq && rocketPair.reverse.saveRecognizedSeq){
				 System.out.println(rocketPair.name +" >>> combining finally recognized seq......");
				 splitedSeqFiles_F=new ArrayList<String>();
				 splitedSeqFiles_R=new ArrayList<String>();
				 for(int i=0;i<splitedRockets12.size();i++){
					splitedSeqFiles_F.add(
					   splitedRockets12.get(i).get(r).forward.recognizedSeqFile
					);
					splitedSeqFiles_R.add(
					   splitedRockets12.get(i).get(r).reverse.recognizedSeqFile					
					);
				 }
				 
				 file=rocketPair.forward.rocketName+".";
				 for(String fileName:splitedSeqFiles_F) {
					 if(fileName!=null && new File(fileName).exists()) {
						file=fileName; 
						break;
					 }
				 }				
				 seqFormat=FileOperation.getFileFormat(file);
				 outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
				 rocketPair.forward.recognizedSeqFile=outName;
				 fileList_forward.add(outName);
				 
				 file2=rocketPair.reverse.rocketName+".";	
				 for(String fileName:splitedSeqFiles_R) {
					 if(fileName!=null && new File(fileName).exists()) {
						file2=fileName; 
						break;
					 }
				 }			
				 seqFormat=FileOperation.getFileFormat(file2);
				 outName2=combinedSeqOut+"/"
				         +file2.substring(file2.lastIndexOf("/")+1,file2.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName2);
				 rocketPair.reverse.recognizedSeqFile=outName2;
				 fileList_reverse.add(outName2);
			 }else {
	    		 if(rocketPair.forward.isRecognized && rocketPair.reverse.isRecognized)	    		   
	    		   System.out.println(rocketPair.name+" >>> Warning: you did not set to save finally recognized file, please check '-finalSave'!!!");
	    		 else
	    		   System.out.println(rocketPair.name+" >>> Warning: No finally recognized seq, or check if your data or configuration is possibly incorrect for it!!!");	  
			 }
		
			 // for recognized Seq trimmed
			 if(rocketPair.forward.isEndTrimmed && rocketPair.reverse.isEndTrimmed){
				 System.out.println(rocketPair.name+" >>> combining finally recognized seq trimmed......");
				 splitedSeqFiles_F=new ArrayList<String>();
				 splitedSeqFiles_R=new ArrayList<String>();
				 for(int i=0;i<splitedRockets12.size();i++){
					splitedSeqFiles_F.add(
					   splitedRockets12.get(i).get(r).forward.recognizedSeqFileTrimmed
					);
					splitedSeqFiles_R.add(
					   splitedRockets12.get(i).get(r).reverse.recognizedSeqFileTrimmed
					);
				 }
				 
				 file=rocketPair.forward.rocketName+".Trimmed.";	
				 for(String fileName:splitedSeqFiles_F) {
					 if(fileName!=null && new File(fileName).exists()) {
						file=fileName; 
						break;
					 }
				 }
				 seqFormat=FileOperation.getFileFormat(file);
				 outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
				 rocketPair.forward.recognizedSeqFileTrimmed=outName;
				 
				 file2=rocketPair.reverse.rocketName+".Trimmed.";
				 for(String fileName:splitedSeqFiles_R) {
					 if(fileName!=null && new File(fileName).exists()) {
						file2=fileName; 
						break;
					 }
				 }
				 seqFormat=FileOperation.getFileFormat(file2);
				 outName2=combinedSeqOut+"/"
				         +file2.substring(file2.lastIndexOf("/")+1,file2.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName2);
				 rocketPair.reverse.recognizedSeqFileTrimmed=outName2;
			 }else {
				 System.out.println(rocketPair.name+" >>> no trimmed file, please the related parameters!!!");
			 }

			 
			 // for recognized Seq masked
			 if(rocketPair.forward.isEndMasked && rocketPair.reverse.isEndMasked){
				 System.out.println(rocketPair.name+" >>> combining finally recognized seq masked......");
				 splitedSeqFiles_F=new ArrayList<String>();
				 splitedSeqFiles_R=new ArrayList<String>();
				 for(int i=0;i<splitedRockets12.size();i++){
					splitedSeqFiles_F.add(
					   splitedRockets12.get(i).get(r).forward.recognizedSeqFileMasked
					);
					splitedSeqFiles_R.add(
					   splitedRockets12.get(i).get(r).reverse.recognizedSeqFileMasked
					);
				 }
				 
				 file=rocketPair.forward.rocketName+".Masked.";	
				 for(String fileName:splitedSeqFiles_F) {
					 if(fileName!=null && new File(fileName).exists()) {
						file=fileName; 
						break;
					 }
				 }
				 seqFormat=FileOperation.getFileFormat(file);
				 outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
				 rocketPair.forward.recognizedSeqFileMasked=outName;
				 
				 file2=rocketPair.reverse.rocketName+".Masked.";	
				 for(String fileName:splitedSeqFiles_R) {
					 if(fileName!=null && new File(fileName).exists()) {
						file2=fileName; 
						break;
					 }
				 }				 	
				 seqFormat=FileOperation.getFileFormat(file2);
				 outName2=combinedSeqOut+"/"
				         +file2.substring(file2.lastIndexOf("/")+1,file2.lastIndexOf("."))
			             +".final."+seqFormat;
				 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName2);
				 rocketPair.reverse.recognizedSeqFileMasked=outName2;
			 }else {
				 System.out.println(rocketPair.name+" >>> no masked file, please the related parameters!!!");
			 }
			 
			 //for component-recognized Seq required to be saved
			 if(!rocketPair.reverse.finalSeqCompo.equalsIgnoreCase(BARCODE_NAME_DEFINITION)){
				for(int c=0; c<rocketPair.reverse.savedCompoAlignedSeqFiles.size();c++) {
				  System.out.println(rocketPair.name+" >>> combining component-recognized seq file required to be saved......");
				  splitedSeqFiles_F=new ArrayList<String>();
				  splitedSeqFiles_R=new ArrayList<String>();
				  for(int i=0;i<splitedRockets12.size();i++){
					splitedSeqFiles_F.add(
					   splitedRockets12.get(i).get(r).forward.savedCompoAlignedSeqFiles.get(c)
					);
					splitedSeqFiles_R.add(
					   splitedRockets12.get(i).get(r).reverse.savedCompoAlignedSeqFiles.get(c)
					);
				  }
				  
				  file=rocketPair.forward.rocketName+".";	
				  for(String fileName:splitedSeqFiles_F) {
					if(fileName!=null && new File(fileName).exists()) {
					   file=fileName; 
					   break;
				    }
				  }
				  seqFormat=FileOperation.getFileFormat(file);
				  outName=combinedSeqOut+"/"
				         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
			             +".required."+seqFormat;
				  SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
				  rocketPair.forward.savedCompoAlignedSeqFiles.set(c,outName);
				 
				  file2=rocketPair.reverse.rocketName+".";	
				  for(String fileName:splitedSeqFiles_R) {
					if(fileName!=null && new File(fileName).exists()) {
					   file2=fileName; 
					   break;
				    }
				  }
				  seqFormat=FileOperation.getFileFormat(file2);
				  outName2=combinedSeqOut+"/"
				         +file2.substring(file2.lastIndexOf("/")+1,file2.lastIndexOf("."))
			             +".required."+seqFormat;
				  SeqOperation.combineSeqFile(splitedSeqFiles_R,outName2);
				  rocketPair.reverse.savedCompoAlignedSeqFiles.set(c,outName2);
				}
			 }
			 
			 seqPairRockets.add(rocketPair);
			 rocketPair=null;
			 splitedSeqFiles_F=null;
			 splitedSeqFiles_R=null;
	     }	
		 
		 String fileListOut_forward=combinedSeqOut+"/fileList_forward.txt";
		 FileOperation.saveList(fileList_forward, null, fileListOut_forward);
		 String fileListOut_reverse=combinedSeqOut+"/fileList_reverse.txt";
		 FileOperation.saveList(fileList_reverse, null, fileListOut_reverse);
		 fileList_forward=null;
		 fileList_reverse=null;

	 }else{
		 System.out.println("Warning: No sequences recognized, please check if your data or configuration is correct!!!");
	 }

	 splitedRockets12=null;
	 saveRecognizedSeq=null;
	 isRecognized=null;
	 isEndTrimmed=null;
	 isEndMasked=null;
	 saveRecognizedSeq2=null;
	 isRecognized2=null;
	 isEndTrimmed2=null;
	 isEndMasked2=null;
	 
	 setSeqPairRockets(seqPairRockets);
	 
	 if(tmpFiles!=null){
		for(String tmpFile: tmpFiles){
		   FileOperation.delFile(tmpFile);
		}
	 }
	 tmpFiles=null;
	 if(forSplitDir!=null) FileOperation.delFolder(forSplitDir);
	 if(revSplitDir!=null) FileOperation.delFolder(revSplitDir);
	 
	 System.out.println("......Seq Recognization Done......");
	 
	 return seqPairRockets;
	 
  } 

  void launchSeqRocket(SeqRocket seqRocket, List<SeqCompoAlignInfo> barNoExactSeqObj, 
		String barNoExactLeftSubSeqFile, String seqOutDir){		
		  	 
	    if(!seqRocket.seqCompoFeatures.compoNames.contains(BARCODE_NAME_DEFINITION)){
		     System.out.println("Warning: "
	           +"You did not set barcode or primer sequence for experiment '"
		       +seqRocket.rocketName+"', so no sequence extraction for this experiment!!!");
			 
		     seqRocket.isDone=false;
		     return;
	    }
	      
		String barBlastCMD="";
		String tarBlastCMD="";
		String baitBlastCMD="";
		String blastCMD="";	
		int blastWordSize=7;
		String blastTask="blastn-short";	
		String tarSeqFile="";
		String tarBaitSeqFile="";
		String querySeqFile="";
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String blastOutFile=tmpDir+"/SeqRocket.BlastOut."+timeStamp+".txt";
		  
		List<SeqCompoAlignInfo> tarSeqObj = null;		
		List<SeqCompoAlignInfo> tarExactSeqObj=null;
		List<SeqCompoAlignInfo> tarNoExactSeqObj=null;
		List<SeqCompoAlignInfo> tarBLASTSeqObj=null;
		List<SeqCompoAlignInfo> combinedTarBLASTSeqObj = null;
		List<SeqCompoAlignInfo> tarNoBLASTSeqObj=null;
		List<ArrayList<SeqCompoAlignInfo>> baitObj=null;
		List<ArrayList<SeqCompoAlignInfo>> rightObj=null;
		List<Integer> blastWordSizeList;
		String seqCompoName="";
		String prevLeftSeqCompoName="";	 
		String prevRightSeqCompoName="";
		List<String> tmpFiles=new ArrayList<String>();
		boolean doBarBLAST=true;		
		int seqNum=0;
		seqNum=SeqOperation.getSeqNum(barNoExactLeftSubSeqFile);		
		barBlastCMD="-evalue 10000 -max_target_seqs "+seqNum;	
		if(seqNum==0) doBarBLAST=false;
		String currName="";
		//String saveName="";
		SeqCompoRecognizer recognizer;
		SeqCompoRecognizer primerCont;
		  
		if(seqRocket.isActive){
			  
		   tarSeqObj=new ArrayList<SeqCompoAlignInfo> ();
		   //currName=seqRocket.rocketName+".";
		   initSeqAlignArray(barNoExactSeqObj,seqRocket.seqCompoFeatures.compoNames.size());						
		   boolean leftTrim=false;
			
		   for(int k=0;k<seqRocket.seqRecognizers.size();k++){	 // for Left seqRecognizer			   
			  recognizer=seqRocket.seqRecognizers.get(k);
			  if(recognizer==null) break;
			  seqCompoName=seqRocket.seqCompoFeatures.compoNames.get(recognizer.index);	  			  
			  //currName=currName+seqCompoName+"-";
			  currName=seqRocket.rocketName+"."+seqCompoName;
			  System.out.println(recognizer.seqName+" are being recognized......");
			   
			  if(recognizer.side.equalsIgnoreCase(SEQ_LEFT_SIDE)){
				   
			    if(seqCompoName.equals(BARCODE_NAME_DEFINITION)){             		  
				    tarExactSeqObj=getSeqObj(recognizer.exactAlignedSeqFile);	 
				    initSeqAlignArray(tarExactSeqObj,seqRocket.seqCompoFeatures.compoNames.size());
				    setSeqLeftExactAlignInfo(tarExactSeqObj,seqRocket,seqCompoName); 
				    tarBLASTSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				    if(doBarBLAST){
					  querySeqFile=recognizer.seqFASTAFile;
					  blastWordSize=recognizer.blastWordSize;
					  blastTask=recognizer.blastTask;
					  blastCMD="blastn -task "+blastTask+" "
					          +barBlastCMD+" -word_size="+blastWordSize
					          +" -query "+querySeqFile
					          +" -subject "+barNoExactLeftSubSeqFile
					          +" -out "+blastOutFile+" -outfmt 6";
					  //System.out.println(blastCMD+"......");
					  SeqOperation.runBLAST(blastCMD);			  
					  tarBLASTSeqObj=getLeftSideBLASTSeq(barNoExactSeqObj,recognizer,
							  seqCompoName,seqRocket,blastOutFile);		  
					  FileOperation.delFile(blastOutFile);				 		  
				    }  
				    prevLeftSeqCompoName=seqCompoName;				  
	                tarSeqObj=combineSeqObj(tarExactSeqObj,tarBLASTSeqObj);
	                System.out.println(seqRocket.rocketName+">>> Got "+seqCompoName
	                     +" Recognized Seq: "+tarSeqObj.size()
				         +" (Exact:"+tarExactSeqObj.size()+", NoExact:"+tarBLASTSeqObj.size()+")");			  
	                tarBLASTSeqObj=null;
			    }else if(seqCompoName.equals(BAIT_NAME_DEFINITION) 
			    		 || seqCompoName.equals(BAITBRK_NAME_DEFINITION) 
					     || seqCompoName.equals(BAITARM_NAME_DEFINITION)){
				    primerCont=null;	
					tarNoBLASTSeqObj=tarSeqObj;  
					combinedTarBLASTSeqObj=new ArrayList<SeqCompoAlignInfo> ();					
	                blastWordSizeList=new ArrayList<Integer> ();
				    if(seqCompoName.equals(BAITARM_NAME_DEFINITION) 
				    		|| seqCompoName.equals(BAITBRK_NAME_DEFINITION)){
				      if(recognizer.seqLength>=2000) blastWordSizeList.add(128);
					  if(recognizer.seqLength>=1000) blastWordSizeList.add(64);
					  if(recognizer.seqLength>=500) blastWordSizeList.add(48);
					  if(recognizer.seqLength>=50) blastWordSizeList.add(21);
				      blastWordSizeList.add(11); 
				    }else if(seqCompoName.equals(BAIT_NAME_DEFINITION)){
					  if(recognizer.seqLength>=500) blastWordSizeList.add(48);
					  if(recognizer.seqLength>=100) blastWordSizeList.add(32);
					  blastWordSizeList.add(21);
					  blastWordSizeList.add(16);
					  blastWordSizeList.add(11);
					  blastWordSizeList.add(7);				  
					  blastWordSizeList.add(4);
					  primerCont=seqRocket.seqRecognizers.get(
						 seqRocket.seqCompoFeatures.compoNames.indexOf(PRIMERCONT_NAME_DEFINITION)
					  );	
				    }
				    
				    for(int c=0;c<blastWordSizeList.size();c++){
					  recognizer.blastWordSize=blastWordSizeList.get(c);
					  recognizer.minAlignLen=recognizer.blastWordSize;
					  if(seqCompoName.equals(BAIT_NAME_DEFINITION) 
							  && blastWordSizeList.get(c)<=primerCont.blastWordSize){
						 recognizer.leftSubForBLAST=true;
						 recognizer.territoryLen=primerCont.territoryLen;
						 recognizer.minAlignLen=primerCont.minAlignLen;
						 recognizer.blastWordSize=primerCont.blastWordSize;
						 recognizer.seqFASTAFile=primerCont.seqFASTAFile;
					  }
					  if(recognizer.leftSubForBLAST && leftTrim){
						tarBaitSeqFile=tmpDir+"/"+seqCompoName+".LeftTrim.LeftSub.ForBLAST."+c+".fna";	   
					  }else if(recognizer.leftSubForBLAST){
						tarBaitSeqFile=tmpDir+"/"+seqCompoName+".LeftSub.ForBLAST."+c+".fna";	   
					  }if(leftTrim){
					    tarBaitSeqFile=tmpDir+"/"+seqCompoName+".LeftTrim.ForBLAST."+c+".fna";	   
					  }else{
						tarBaitSeqFile=tmpDir+"/"+seqCompoName+".ForBLAST."+c+".fna";	   
					  }
					  createLeftRecognizerSeq(tarNoBLASTSeqObj,seqRocket,recognizer,tarBaitSeqFile);
	                  tmpFiles.add(tarBaitSeqFile);				  
					  blastWordSize=recognizer.blastWordSize;
					  if(blastWordSize<=16 || recognizer.seqLength<50){	  
						recognizer.blastTask="blastn-short";
					  }else if(blastWordSize<=24){
						recognizer.blastTask="blastn";
					  }else{
					    recognizer.blastTask="megablast";
					  }
					  blastTask=recognizer.blastTask;
					  baitBlastCMD="blastn -task "+blastTask+" -word_size="+blastWordSize+" -evalue 10000";
					  seqNum=SeqOperation.getSeqNum(tarBaitSeqFile);		
					  baitBlastCMD=baitBlastCMD+" -max_target_seqs "+seqNum;         	  
					  querySeqFile=recognizer.seqFASTAFile;
					  blastCMD=baitBlastCMD+" -query "+querySeqFile+" -subject "
					          +tarBaitSeqFile+" -out "+blastOutFile+" -outfmt 6";						 
					  if(seqNum>0){
						//System.out.println(blastCMD+"......");
						SeqOperation.runBLAST(blastCMD);		
						baitObj=getBaitBLASTSeq(tarNoBLASTSeqObj,recognizer,seqCompoName,seqRocket,blastOutFile);
						FileOperation.delFile(blastOutFile);
						tarBLASTSeqObj=baitObj.get(0);   //seqObj having BLAST alignment,meaning those seq contain bait, breaksite or bait arm region.
						tarNoBLASTSeqObj=baitObj.get(1); //seqObj not having BLAST alignment, meaning those seq don't contain bait, breaksite or bait arm region.
						if(tarBLASTSeqObj.size()>0) 
						   combinedTarBLASTSeqObj=combineSeqObj(combinedTarBLASTSeqObj,tarBLASTSeqObj);						
					  }else{
						//System.out.println(" No Sequence for "+recognizer.seqName+" Blast!");
						break;
					  }	  
					  baitObj=null;				  	
				    }// for blastWordSize
				    
				    tarSeqObj=combinedTarBLASTSeqObj;
				    if(tarNoBLASTSeqObj.size()>0){			    
					  //Set alignArray info for this recognizer for those seqObj without BLAST alignment by copying the align position of previous recognizer in order to show that this seq junction goes to out of bait region.
				      copyAlignInfo(tarNoBLASTSeqObj,seqRocket,prevLeftSeqCompoName,seqCompoName);
					  tarSeqObj=combineSeqObj(combinedTarBLASTSeqObj,tarNoBLASTSeqObj);
				    }
				    
				    System.out.println(seqRocket.rocketName+">>>Seq harbouring portion of "+seqCompoName
					   +" / totally recognized seq: "+combinedTarBLASTSeqObj.size() +"/"+tarSeqObj.size());
				    
				    tarBLASTSeqObj=null;
				    combinedTarBLASTSeqObj=null;
				    tarNoBLASTSeqObj=null;
					   
				    prevLeftSeqCompoName=seqCompoName;			  
			  
			    }else{
				    setSeqLeftExactAlignInfo(tarSeqObj,seqRocket,seqCompoName); 		  
				    tarExactSeqObj=getRecognizedSeq(tarSeqObj,seqRocket,seqCompoName);
				    tarNoExactSeqObj=getNoRecognizedSeq(tarSeqObj,seqRocket,seqCompoName);
				    tarSeqObj=null;		
				    tarBLASTSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				    if(recognizer.leftSubForBLAST && leftTrim){
					  tarSeqFile=tmpDir+"/"+seqCompoName+".Trim.leftsub.ForBLAST.fna";	   
				    }else if(recognizer.leftSubForBLAST){
					  tarSeqFile=tmpDir+"/"+seqCompoName+".LeftSub.ForBLAST.fna";	   
				    }else if(leftTrim){
				      tarSeqFile=tmpDir+"/"+seqCompoName+".Trim.ForBLAST.fna";	
				    }else{
					  tarSeqFile=tmpDir+"/"+seqCompoName+".ForBLAST.fna";	    
				    }	
				    createLeftRecognizerSeq(tarNoExactSeqObj,seqRocket,recognizer,tarSeqFile);	
	                tmpFiles.add(tarSeqFile);		  
				    blastWordSize=recognizer.blastWordSize;
	                blastTask=recognizer.blastTask;			  
				    tarBlastCMD="blastn -task "+blastTask+" -word_size="+blastWordSize+" -evalue 10000";
				    seqNum=SeqOperation.getSeqNum(tarSeqFile);		
				    tarBlastCMD=tarBlastCMD+" -max_target_seqs "+seqNum;         	  
				    querySeqFile=recognizer.seqFASTAFile;
				    blastCMD=tarBlastCMD+" -query "+querySeqFile
						  +" -subject "+tarSeqFile+" -out "+blastOutFile+" -outfmt 6";				 
				    if(seqNum>0){
					  //System.out.println(blastCMD+"......");
				      SeqOperation.runBLAST(blastCMD);		
				      tarBLASTSeqObj=getLeftSideBLASTSeq(
				    		tarNoExactSeqObj,recognizer,seqCompoName,seqRocket,blastOutFile
				      );          
					  FileOperation.delFile(blastOutFile);			
				    }else{
					//System.out.println("No Sequence for "+recognizer.seqName+" Blast!");
				    }	 
				    tarNoExactSeqObj=null;	 
				    tarSeqObj=combineSeqObj(tarExactSeqObj,tarBLASTSeqObj);
				    prevLeftSeqCompoName=seqCompoName;
				    System.out.println(seqRocket.rocketName+" >>> Got "+seqCompoName
					  +" Recognized Seq: "+tarSeqObj.size()
					  +" (Exact: "+tarExactSeqObj.size()+", NoExact: "+tarBLASTSeqObj.size()+")");				  
			    }
			  
			    tarExactSeqObj=null;
			    tarBLASTSeqObj=null;
			     
			    seqRocket.seqRecognizers.get(k).done=true;						  
			    if(recognizer.leftShiftForNext) leftTrim=true;				    

			    if(seqRocket.saveRecognizedSeqAsHTML)
			       recognizer.alignedSeqHTMLFile=seqOutDir+currName+".AlignedSeq.html";
			    
		    	recognizer.alignedSeqFileMasked=seqOutDir+currName+".AlignedSeqMasked.fna";
		    	recognizer.alignedSeqFileTrimmed=seqOutDir+currName+".AlignedSeqTrimmed.fna";	
		    	recognizer.alignedSeqFile=seqOutDir+currName+".AlignedSeq";		    	
						  
			    if(recognizer.saveAsFASTQ){
			   	   recognizer.alignedSeqFile=recognizer.alignedSeqFile+".fastq";
				   saveAsFASTQFile(tarSeqObj,recognizer.alignedSeqFile);	
				   recognizer.isFASTQSaved=true;					
				   seqRocket.savedCompoAlignedSeqFiles.add(recognizer.alignedSeqFile);
				   //seqRocket.recognizedSeqFile=recognizer.alignedSeqFile;
				   //seqRocket.isRecognized=true;				  				    
				}else if(recognizer.saveAsFASTA){
				   recognizer.alignedSeqFile=recognizer.alignedSeqFile+".fna";
				   saveAsFASTAFile(tarSeqObj,recognizer.alignedSeqFile);	
				   recognizer.isFASTASaved=true;										  
				   seqRocket.savedCompoAlignedSeqFiles.add(recognizer.alignedSeqFile);
				   //seqRocket.recognizedSeqFile=recognizer.alignedSeqFile;
				   //seqRocket.isRecognized=true;			
				}
				
			    seqRocket.finalSeqCompo=seqCompoName;
				seqRocket.isDone=true;				   
				System.out.println(seqRocket.rocketName +" >>> Recognized "+currName+" is saved!!!");				
			   
			  } // if READ_LEFT_SIDE
			  recognizer=null;	  
		    } // for Left seqRecognizer
			
		    if(doSeqRight){
		      prevRightSeqCompoName=null;
			  boolean rightTrim=false;
			  for(int k=seqRocket.seqRecognizers.size()-1;k>=0;k--){	 // for Right seqRecognizer	
				 recognizer=seqRocket.seqRecognizers.get(k);
				 if(recognizer==null) break;
				 seqCompoName=seqRocket.seqCompoFeatures.compoNames.get(recognizer.index);
				 if(recognizer.side.equalsIgnoreCase(SEQ_RIGHT_SIDE)){
					tarExactSeqObj=new ArrayList<SeqCompoAlignInfo> ();
					tarNoExactSeqObj=new ArrayList<SeqCompoAlignInfo> ();	
					setSeqRightExactAlignInfo(tarSeqObj,seqRocket,seqCompoName); 		  
					tarExactSeqObj=getRecognizedSeq(tarSeqObj,seqRocket,seqCompoName);
					tarNoExactSeqObj=getNoRecognizedSeq(tarSeqObj,seqRocket,seqCompoName);
					tarNoBLASTSeqObj=tarNoExactSeqObj;
					combinedTarBLASTSeqObj=new ArrayList<SeqCompoAlignInfo> ();
					tarSeqObj=null;							  
						  
			        blastWordSizeList=new ArrayList<Integer> ();                   
					if(recognizer.seqLength>=500) blastWordSizeList.add(48);
					if(recognizer.seqLength>=100) blastWordSizeList.add(32);
					if(recognizer.seqLength>=50) blastWordSizeList.add(21);
					blastWordSizeList.add(16);
					blastWordSizeList.add(11);
					blastWordSizeList.add(7);					
											    
					for(int c=0;c<blastWordSizeList.size();c++){
					   recognizer.blastWordSize=blastWordSizeList.get(c);
					   //recognizer.minAlignLen=recognizer.blastWordSize;
					   if(recognizer.rightSubForBLAST && rightTrim){
						  tarSeqFile=tmpDir+"/"+seqCompoName+".Trim.RightSub.ForBLAST."+c+".fna";	   
					   }else if(recognizer.rightSubForBLAST){
						  tarSeqFile=tmpDir+"/"+seqCompoName+".RightSub.ForBLAST."+c+".fna";	     
					   }else if(rightTrim){
						  tarSeqFile=tmpDir+"/"+seqCompoName+".Trim.ForBLAST."+c+".fna";	   
					   }else{
						  tarSeqFile=tmpDir+"/"+seqCompoName+".ForBLAST."+c+".fna";
					   }	
					   createRightRecognizerSeq(tarNoBLASTSeqObj,seqRocket,recognizer,tarSeqFile);	
			           tmpFiles.add(tarSeqFile);				  
					   blastWordSize=recognizer.blastWordSize;
					   if(blastWordSize<=16 || recognizer.seqLength<50){	  
						  recognizer.blastTask="blastn-short";
					   }else if(blastWordSize<=24){
						  recognizer.blastTask="blastn";
					   }else{
						  recognizer.blastTask="megablast";
					   }
					   blastTask=recognizer.blastTask;
					   tarBlastCMD="blastn -task "+blastTask+" -word_size="+blastWordSize+" -evalue 10000";
					   seqNum=SeqOperation.getSeqNum(tarSeqFile);		
					   tarBlastCMD=tarBlastCMD+" -max_target_seqs "+seqNum;         	  
					   querySeqFile=recognizer.seqFASTAFile;
				       //blastOutFile=blastOutFile.substring(0,blastOutFile.lastIndexOf("."))+"."+c+".out";
					   blastCMD=tarBlastCMD+" -query "+querySeqFile+" -subject "
							          +tarSeqFile+" -out "+blastOutFile+" -outfmt 6";	
					   rightObj=new ArrayList<ArrayList<SeqCompoAlignInfo>>();
					   if(seqNum>0){
						  //System.out.println(blastCMD+"......");
						  SeqOperation.runBLAST(blastCMD);		
						  rightObj=getRightSideBLASTSeq(tarNoBLASTSeqObj,recognizer,seqCompoName,seqRocket,blastOutFile);
						  FileOperation.delFile(blastOutFile);
						  tarBLASTSeqObj=rightObj.get(0);   //seqObj having BLAST alignment,meaning those seq contain or partially contain read right recognizer.
						  tarNoBLASTSeqObj=rightObj.get(1); //seqObj not having BLAST alignment, meaning those seq don't contain read right recognizer.
						  if(tarBLASTSeqObj.size()>0) 
							 combinedTarBLASTSeqObj=combineSeqObj(combinedTarBLASTSeqObj,tarBLASTSeqObj);						
					   }else{
						 //System.out.println(" No Sequence for "+recognizer.seqName+" Blast!");
						  break;
					   }	  
					   rightObj=null;				  	
					}// for blastWordSize				    
				    
					tarNoExactSeqObj=combinedTarBLASTSeqObj;
					if(tarNoBLASTSeqObj.size()>0){			    
					   //Set alignArray info of this recognizer for those seqObj without BLAST alignment by shifting align position of previous recognizer.
					   copyAlignInfo(tarNoBLASTSeqObj,seqRocket,prevRightSeqCompoName,seqCompoName);
					   tarNoExactSeqObj=combineSeqObj(combinedTarBLASTSeqObj,tarNoBLASTSeqObj);
					}
				 
					tarSeqObj=combineSeqObj(tarExactSeqObj,tarNoExactSeqObj);
			          
					System.out.println(seqRocket.rocketName+" >>> Seq harbouring portion of "+seqCompoName
							+" / totally recognized seq: "
							+(tarExactSeqObj.size()+combinedTarBLASTSeqObj.size()) +"/"+tarSeqObj.size());
					   
					tarBLASTSeqObj=null;
					combinedTarBLASTSeqObj=null;					 
					tarNoBLASTSeqObj=null;
					tarNoExactSeqObj=null;
					tarExactSeqObj=null;
					   
					prevRightSeqCompoName=seqCompoName;	
					 
					seqRocket.seqRecognizers.get(k).done=true;						  
					if(recognizer.rightShiftForNext) rightTrim=true;	

			     } // if READ_RIGHT_SIDE
	
			  }// for Right seqRecognizer
			} // do seq right side		    
		 
		    seqRocket.isRecognized=true;
		    if(seqRocket.saveRecognizedSeq) saveRecognizedSeq(tarSeqObj, seqRocket);
		    if(seqRocket.saveRecognizedSeqAsHTML) saveRecognizedSeqAsHTML(tarSeqObj, seqRocket);
		    
			tarSeqObj=null; //release memory space			

		}// is rocket active?  
		
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
		   FileOperation.delFile(tmpFile);
		}
	    tmpFiles=null;

   }
 
}
  
  
