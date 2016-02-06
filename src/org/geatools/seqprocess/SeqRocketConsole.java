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
import java.io.*;  
import java.text.SimpleDateFormat;
import java.util.ArrayList; 
import java.util.Calendar;
import java.util.List;

import org.geatools.data.structure.BlastInfo;
import org.geatools.data.structure.SeqCompoAlignInfo;
import org.geatools.data.structure.SeqCompo;
import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqCompoType;
import org.geatools.data.structure.SeqRocket;
import org.geatools.data.structure.SeqRocketPair;
import org.geatools.operation.FileOperate;

public class  SeqRocketConsole {
  
  List<SeqRocket> seqRockets;
  List<SeqRocketPair> seqPairRockets;
  boolean isRocketsOK=false;
  
  BlastInfo blastInfoObj=new BlastInfo();  
  double extendTerritoryRatio=0.3d;
  
  int barTerritoryLen=12;
  int minBarcodeLen=8;
  int minRedPrimerLen=minBarcodeLen;
  
  double barMinAlignRatio=0.7d;
  double barMaxMismatchRatio=0.1d;
  double barMaxGapRatio=0.1d;
 
  double primerMinAlignRatio=0.7d;
  double primerMaxMismatchRatio=0.1d;
  double primerMaxGapRatio=0.1d;
  int primerTrimLeftShift=5;
  
  double primerContMinAlignRatio=0.8d;
  double primerContMaxMismatchRatio=0.1d;
  double primerContMaxGapRatio=0.1d;
  
  double baitMinAlignRatio=0.7d;
  double baitMaxMismatchRatio=0.25d;
  double baitMaxGapRatio=0.25d;
  int baitTrimLeftShift=10;
  
  int baitBrkMinAlign=12;
  double baitBrkMinAlignRatio=0.7d;
  double baitBrkMaxMismatchRatio=0.25d;
  double baitBrkMaxGapRatio=0.25d;
  int baitBrkTrimLeftShift=10;
  
  int baitArmMinAlign=20;
  double baitArmMinAlignRatio=0.7d;
  double baitArmMaxMismatchRatio=0.25d;
  double baitArmMaxGapRatio=0.25d;
  int baitArmTrimLeftShift=10;

  //String resultDir="";
 
  String barcodeName="Barcode";
  String primerName="Primer";
  String primerContName="PrimerCont";
  String baitName="Bait";
  String baitBrkName="BaitBrk";
  String baitArmName="BaitArm";
  String bodyName="Body";
  String rPrimerName="ReversePrimer";
  
  String barcodeColor="black";
  String primerColor="red";
  String primerContColor="yellow";
  String baitColor="pink";
  String baitBrkColor="gray";
  String baitArmColor="brown";
  String bodyColor="green";
  String rPrimerColor="blue";
  
  
  static String homeDir=".";
  List<String> tmpFiles;
  String tmpDir;
  String dataDir;
  
  String regexSeqID = "(\\s+SeqID=)(\\d+)(#\\s*)";
  
  public SeqRocketConsole(){ 

  }
  public SeqRocketConsole(String homeDir,String dataDir,String tmpDir){ 
	  setHomeDir(homeDir);
	  setDataDir(dataDir);
	  setTmpDir(tmpDir);
  } 
  
  public void setHomeDir(String dir){
      homeDir=dir;
  }
  
  public void setDataDir(String dir){
	  dataDir=dir;
  }
  
  public void setTmpDir(String dir){
	  tmpDir=dir;
  }
  
  public boolean isSeqRocketsOK(){
	  return isRocketsOK;
  }
  public void setSeqRockets(List<SeqRocket> rockets){
	 seqRockets=rockets;
	 if(rockets!=null && rockets.size()>0) isRocketsOK=true;
	 else isRocketsOK=false;;
  }
  public List<SeqRocket>  getSeqRockets(){
	 return seqRockets;
  }  
  public void setSeqPairRockets(List<SeqRocketPair> pairRockets){
	 seqPairRockets=pairRockets;
	 if(pairRockets!=null && pairRockets.size()>0) isRocketsOK=true;
	 else isRocketsOK=false;;
  }
  public List<SeqRocketPair>  getSeqPairRockets(){
	 return seqPairRockets;
  }
 
  void configRecognizer( List<SeqRocket> rockets){
    
	SeqCompoRecognizer barcode;
	SeqCompoRecognizer primer;
	SeqCompoRecognizer primerCont;
	SeqCompoRecognizer bait;
	SeqCompoRecognizer baitBrk;
	SeqCompoRecognizer baitArm;
	int barRawSeqLen=0;
	for(int i=0;i<rockets.size();i++){
	  
      barcode=null;
      primer=null;
      primerCont=null;
	  bait=null;
      baitBrk=null;	
      baitArm=null;
	  
      SeqCompoRecognizer recognizer;
	  String seqType="";
	  for(int k=0;k<rockets.get(i).seqRecognizer.size();k++){
	   recognizer=rockets.get(i).seqRecognizer.get(k);	
	   seqType=rockets.get(i).seqTypeInfo.seqTypeName.get(recognizer.index);
	   if(seqType.equals(barcodeName))
		 barcode=recognizer;
	   if(seqType.equals(primerName))
		 primer=recognizer;
	   if(seqType.equals(primerContName))
		 primerCont=recognizer;
	   if(seqType.equals(baitName))
	     bait=recognizer;
	   if(seqType.equals(baitBrkName))
	     baitBrk=recognizer;
	   if(seqType.equals(baitArmName))
	     baitArm=recognizer;
	  }
	  recognizer=null;
	  
	  if(barcode!=null){

		  if(barcode.rawSeq==null){
			barcode.seq=primer.rawSeq.substring(0,minBarcodeLen);
			barRawSeqLen=0;
		  }else if(barcode.rawSeq.length()<minBarcodeLen){
			barcode.seq=barcode.rawSeq+primer.rawSeq.substring(0,minBarcodeLen-barcode.rawSeq.length());
			barRawSeqLen=barcode.rawSeq.length();
		  }else{
			barcode.seq=barcode.rawSeq;
			barRawSeqLen=barcode.rawSeq.length();
		  }
		  barcode.seqLength=barcode.seq.length();
		  barcode.territoryLen=barcode.seqLength+(int) Math.ceil(barcode.seqLength*extendTerritoryRatio);	
		  barcode.maxExactStart=1;	
		  barcode.minAlignRatio=barMinAlignRatio;
          barcode.maxMismatchRatio=barMaxMismatchRatio;
          barcode.maxGapRatio=barMaxGapRatio;			
		  barcode.minAlignLen=(int) Math.ceil(barcode.seqLength*barcode.minAlignRatio);
		  barcode.maxMismatchNum=(int) Math.ceil(barcode.seqLength*barcode.maxMismatchRatio);
		  barcode.maxGapNum=(int) Math.ceil(barcode.seqLength*barcode.maxGapRatio);	
		  if(barcode.minAlignLen<=12){	  
			barcode.blastWordSize=4;
			barcode.blastTask="blastn-short";
		  }else if(barcode.minAlignLen<=21){
			barcode.blastWordSize=7;
			barcode.blastTask="blastn-short";
		  }else if(barcode.minAlignLen<=50){
			barcode.blastWordSize=11;
			barcode.blastTask="blastn-short";
	      }else if(barcode.minAlignLen<=75){
			barcode.blastWordSize=16;
			barcode.blastTask="blastn";
	      }else{
		    barcode.blastWordSize=24;
			barcode.blastTask="megablast";
		  }
		  
		  if(barRawSeqLen==0){
			barcode.maxQStart=barcode.seqLength/2;
			barcode.maxSStart=1;
		  }else if(barRawSeqLen<4){
			barcode.maxQStart=1;
			barcode.maxSStart=1;
		  }else if(barRawSeqLen<minBarcodeLen){
			barcode.maxQStart=barRawSeqLen/2;
			barcode.maxSStart=barRawSeqLen/2;
          }else{			
			barcode.maxQStart=barcode.seqLength-barcode.minAlignLen+1;
			barcode.maxSStart=barcode.territoryLen-barcode.minAlignLen+1;
		  }	   
		  barcode.leftSubForBlast=false;
		  barcode.leftShiftForNext=false;
		  barcode.saveRecognizedSeq=true;
		  barcode.saveSeqAsFASTAFormat=true;		
		  barcode.leftTrimSave=false;
		  barcode.rightTrimSave=false;
		  barcode.exactAlignedSeqFile=tmpDir+"/"+barcode.seqName+"_ExactAlignedSeq";
		  barcode.tagSeqFile=tmpDir+"/recognizer/barcode/"+barcode.seqName+".fna";	
		  createFASTASeq(barcode.seq,barcode.seqName,barcode.tagSeqFile); 

          rockets.get(i).seqRecognizer.set(barcode.index,barcode);	
          
		  tmpFiles.add(barcode.exactAlignedSeqFile);
		  tmpFiles.add(barcode.tagSeqFile);
      }
      
	  if(primer!=null){		
		  primer.seq=primer.rawSeq;			
		  primer.seqLength=primer.seq.length();
		  if(barcode.leftShiftForNext)
			primer.territoryLen=primer.seqLength;
		  else
			primer.territoryLen=barRawSeqLen+primer.seqLength;
		   
		  primer.territoryLen=primer.territoryLen+(int) Math.ceil(primer.territoryLen*extendTerritoryRatio);	
		  primer.maxExactStart=primer.territoryLen-primer.seqLength+1; 
		  primer.minAlignRatio=primerMinAlignRatio;	
		  primer.maxMismatchRatio=primerMaxMismatchRatio;
		  primer.maxGapRatio=primerMaxGapRatio;		  
		  primer.minAlignLen=(int) Math.ceil(primer.seqLength*primer.minAlignRatio);
		  primer.maxQStart=primer.seqLength-primer.minAlignLen+1;
		  primer.maxSStart=primer.territoryLen-primer.minAlignLen+1;	
		  if(primer.minAlignLen<=12){	  
			  primer.blastWordSize=4;
			  primer.blastTask="blastn-short";
		  }else if(primer.minAlignLen<=21){
			  primer.blastWordSize=7;
			  primer.blastTask="blastn-short";
		  }else if(primer.minAlignLen<=50){
			  primer.blastWordSize=11;
			  primer.blastTask="blastn-short";
		  }else if(primer.minAlignLen<=75){
			  primer.blastWordSize=16;
			  primer.blastTask="blastn";
		  }else{
			  primer.blastWordSize=24;	
			  primer.blastTask="megablast";
		  }	
		  
		  primer.leftSubForBlast=true;
		  primer.leftShiftForNext=true;
		  if(primerCont!=null){
			  primer.trimLeftShift=primerTrimLeftShift;
			  primer.saveRecognizedSeq=false;
			  primer.alignHTMLSave=false;
			  primer.leftTrimSave=false;
			  primer.rightTrimSave=false;
		  }else{
			  primer.trimLeftShift=0;
			  primer.saveRecognizedSeq=true;
			  primer.saveSeqAsFASTAFormat=true;
			  primer.alignHTMLSave=false;
			  if(bait==null)  primer.alignHTMLSave=true;			 
			  primer.leftMaskSave=true;	
			  primer.leftTrimSave=false;
			  primer.rightTrimSave=false;
		  }
		  primer.tagSeqFile=tmpDir+"/recognizer/primer/"+primer.seqName+".fna";
		  createFASTASeq(primer.seq,primer.seqName,primer.tagSeqFile);
		  
		  rockets.get(i).seqRecognizer.set(primer.index,primer);
		  
		  tmpFiles.add(primer.tagSeqFile);
		  
      }	
	  
	  if(primerCont!=null){	
          if(primerCont.rawSeq==null || primerCont.rawSeq.equals("")){
        	  primerCont.rawSeq=bait.rawSeq.substring(0,primer.trimLeftShift);		  
    	  }   
          if(primer.leftShiftForNext){		  
        	  primerCont.seq=primer.seq.substring(
		    		           primer.seqLength-primer.trimLeftShift,primer.seqLength
		                     )+primerCont.rawSeq;	
          }else{
        	  primerCont.seq=primer.seq+primerCont.rawSeq;
		  }
		  
          primerCont.seqLength=primerCont.seq.length();
          primerCont.territoryLen=primerCont.seqLength;
          primerCont.territoryLen=primerCont.territoryLen
        		  +(int) Math.ceil(primerCont.territoryLen*extendTerritoryRatio);	
          primerCont.maxExactStart=primerCont.territoryLen-primerCont.seqLength+1;
          primerCont.minAlignRatio=primerContMinAlignRatio;
          primerCont.maxMismatchRatio=primerContMaxMismatchRatio;
          primerCont.maxGapRatio=primerContMaxGapRatio;
          primerCont.minAlignLen=(int) Math.ceil(primerCont.seqLength*primerCont.minAlignRatio);
		 		 
		  if(primerCont.minAlignLen<=12){	  
			  primerCont.blastWordSize=4;
			  primerCont.blastTask="blastn-short";
		  }else if(primerCont.minAlignLen<=21){
			  primerCont.blastWordSize=7;
			  primerCont.blastTask="blastn-short";
		  }else if(primerCont.minAlignLen<=50){
			  primerCont.blastWordSize=11;
			  primerCont.blastTask="blastn-short";
		  }else if(primerCont.minAlignLen<=75){
			  primerCont.blastWordSize=16;
			  primerCont.blastTask="blastn";
		  }else{
			  primerCont.blastWordSize=24;	
			  primerCont.blastTask="megablast";
		  }	
		 
		  primerCont.maxQStart=primerCont.seqLength-primerCont.minAlignLen+1;
		  primerCont.maxSStart=primerCont.territoryLen-primerCont.minAlignLen+1;	 		  
		  primerCont.leftSubForBlast=true;
		  primerCont.leftShiftForNext=false;
		  if(bait!=null){
			  primerCont.trimLeftShift=primerCont.rawSeq.length();
			  primerCont.saveRecognizedSeq=false;
			  primerCont.alignHTMLSave=false;
			  primerCont.leftTrimSave=false;
			  primerCont.rightTrimSave=false;
          }else{
        	  primerCont.trimLeftShift=0;
        	  primerCont.saveRecognizedSeq=true;
        	  primerCont.saveSeqAsFASTAFormat=true;
        	  primerCont.alignHTMLSave=true;        
        	  primerCont.leftMaskSave=true;
        	  primerCont.leftTrimSave=false;
        	  primerCont.rightTrimSave=false;
              	  
		  }		  
		  primerCont.tagSeqFile=tmpDir+"/recognizer/primer/"+primerCont.seqName+".fna";
		  createFASTASeq(primerCont.seq,primerCont.seqName,primerCont.tagSeqFile);
		 
		  rockets.get(i).seqRecognizer.set(primerCont.index,primerCont);
		  
		  tmpFiles.add(primerCont.tagSeqFile);
		  
      }	

      if(bait!=null){		
		  bait.seq=bait.rawSeq;			  
		  if(primer.leftShiftForNext){	 
	        bait.seq=primer.seq.substring(
	        		  primer.seqLength-primer.trimLeftShift,primer.seqLength
	        		)+bait.seq;	
	      }else if(barcode.leftShiftForNext){
	        bait.seq=primer.seq+bait.seq;
	      }else{
	        bait.seq=barcode.rawSeq+primer.seq+bait.seq;
	      }
	      bait.seqLength=bait.seq.length();
          bait.territoryLen=bait.seqLength+(int) Math.ceil(bait.seqLength*extendTerritoryRatio);	
		  bait.maxExactStart=bait.territoryLen-bait.seqLength+1;		  
          bait.minAlignRatio=baitMinAlignRatio;		  
		  bait.maxMismatchRatio=baitMaxMismatchRatio;
		  bait.maxGapRatio=baitMaxGapRatio;
		  bait.minAlignLen=primerCont.minAlignLen;
		  bait.maxQStart=primerCont.maxQStart;
	      bait.maxSStart=primerCont.maxSStart;	
          bait.blastWordSize=primerCont.blastWordSize;		  
	      
	      bait.leftSubForBlast=false;
		  if(bait.territoryLen<100)  bait.leftSubForBlast=true;
		  
	      bait.leftShiftForNext=true;		  
		  if(baitBrk!=null || baitArm!=null ){
		    bait.trimLeftShift=0;
		    bait.saveRecognizedSeq=false;
		    bait.alignHTMLSave=false;
	        bait.leftTrimSave=false;
	        bait.rightTrimSave=false;
		    bait.leftMaskSave=false;		  
		  }else{
		    bait.trimLeftShift=baitTrimLeftShift;
			bait.saveRecognizedSeq=true;
			bait.saveSeqAsFASTAFormat=true;
			bait.alignHTMLSave=false;					
			bait.leftMaskSave=true;
			bait.leftTrimSave=false;	
			bait.rightTrimSave=false;
		  }
	      bait.tagSeqFile=tmpDir+"/recognizer/bait/"+bait.seqName+".fna";	
	      createFASTASeq(bait.seq,bait.seqName,bait.tagSeqFile);
		  
		  rockets.get(i).seqRecognizer.set(bait.index,bait);
		  
		  tmpFiles.add(bait.tagSeqFile);
		 
      }
	  
	  if(baitBrk!=null){		
		  baitBrk.seq=baitBrk.rawSeq;
	      baitBrk.seqLength=baitBrk.seq.length();
          baitBrk.territoryLen=baitBrk.seqLength+(int) Math.ceil(baitBrk.seqLength*extendTerritoryRatio);	
		  baitBrk.maxExactStart=baitBrk.territoryLen-baitBrk.seqLength+1;
          baitBrk.minAlignRatio=baitBrkMinAlignRatio;		  
		  baitBrk.maxMismatchRatio=baitBrkMaxMismatchRatio;
		  baitBrk.maxGapRatio=baitBrkMaxGapRatio;
		  baitBrk.minAlignLen=baitBrkMinAlign;		
		  baitBrk.maxQStart=baitBrk.seqLength-baitBrk.minAlignLen+1;
	      baitBrk.maxSStart=baitBrk.territoryLen-baitBrk.minAlignLen+1;
		  if(baitBrk.minAlignLen<=12)	  
			baitBrk.blastWordSize=4;
		  else if(baitBrk.minAlignLen<=21)
			baitBrk.blastWordSize=7;
		  else if(baitBrk.minAlignLen<=33)
			baitBrk.blastWordSize=11;
		  else
		    baitBrk.blastWordSize=16;
		
		  baitBrk.leftSubForBlast=false;
		  if(baitBrk.territoryLen<100)  baitBrk.leftSubForBlast=true;
		  
		  baitBrk.leftShiftForNext=true;
		  if(baitArm!=null){
			  baitBrk.trimLeftShift=0;
			  baitBrk.saveRecognizedSeq=false;
			  baitBrk.alignHTMLSave=false;
			  baitBrk.leftTrimSave=false;
			  baitBrk.rightTrimSave=false;
			  baitBrk.leftMaskSave=false;
		  }else{
			  baitBrk.trimLeftShift=baitBrkTrimLeftShift;
			  baitBrk.saveRecognizedSeq=true;
			  baitBrk.saveSeqAsFASTAFormat=true;
			  baitBrk.alignHTMLSave=true;
			  baitBrk.leftMaskSave=true;
			  baitBrk.leftTrimSave=false;
			  baitBrk.rightTrimSave=false;
          }		  
	      baitBrk.tagSeqFile=tmpDir+"/recognizer/bait/"+baitBrk.seqName+".fna";	
	      createFASTASeq(baitBrk.seq,baitBrk.seqName,baitBrk.tagSeqFile);
		  
		  rockets.get(i).seqRecognizer.set(baitBrk.index,baitBrk);
		  
		  tmpFiles.add(baitBrk.tagSeqFile);
		 
      }
	  
	  if(baitArm!=null){		
		  baitArm.seq=baitArm.rawSeq;
	      baitArm.seqLength=baitArm.seq.length();
          baitArm.territoryLen=baitArm.seqLength+(int) Math.ceil(baitArm.seqLength*extendTerritoryRatio);	
		  baitArm.maxExactStart=baitArm.territoryLen-baitArm.seqLength+1;
          baitArm.minAlignRatio=baitArmMinAlignRatio;		  
		  baitArm.maxMismatchRatio=baitArmMaxMismatchRatio;
		  baitArm.maxGapRatio=baitArmMaxGapRatio;
		  baitArm.minAlignLen=baitArmMinAlign;			
		  baitArm.maxQStart=baitArm.seqLength-baitArm.minAlignLen+1;
	      baitArm.maxSStart=baitArm.territoryLen-baitArm.minAlignLen+1;	
          if(baitArm.minAlignLen<=12)	  
			baitArm.blastWordSize=4;
		  else if(baitArm.minAlignLen<=21)
			baitArm.blastWordSize=7;
		  else if(baitArm.minAlignLen<=33)
			baitArm.blastWordSize=11;
		  else
		    baitArm.blastWordSize=16;		
			
	      baitArm.leftSubForBlast=false;
		  if(baitArm.territoryLen<100)  baitArm.leftSubForBlast=true;
		  
	      baitArm.leftSubForBlast=true;
	      baitArm.leftShiftForNext=true;
		  baitArm.trimLeftShift=baitArmTrimLeftShift;
	      baitArm.saveRecognizedSeq=true;
		  baitArm.saveSeqAsFASTAFormat=true;
		  baitArm.alignHTMLSave=true;
		  baitArm.leftMaskSave=true;
	      baitArm.leftTrimSave=false;
	      baitArm.rightTrimSave=false;
		
	      baitArm.tagSeqFile=tmpDir+"/recognizer/bait/"+baitArm.seqName+".fna";	
	      createFASTASeq(baitArm.seq,baitArm.seqName,baitArm.tagSeqFile);
		  
		  rockets.get(i).seqRecognizer.set(baitArm.index,baitArm);
		  
		  tmpFiles.add(baitArm.tagSeqFile);
		 
      }
	 	  
	  barcode=null;
	  primer=null;
	  primerCont=null;
	  bait=null; 
	  baitBrk=null;
	  baitArm=null;
	
	}	
   
  }
  
  void createFASTASeq(String seq, String seqName, String outSeqFile){
  
	try{    
        BufferedWriter writer=null;
        writer=new BufferedWriter(new FileWriter(outSeqFile));
		String nameLine;
	    //int seqID=0;
	
        nameLine=">"+seqName+" "+seq.length();
		writer.write(nameLine);
		writer.newLine();
		  //writer.flush(); 
		writer.write(seq);
		writer.newLine();
		writer.flush(); 
		
		writer.close();
		
	}catch(IOException e){
        System.out.println(e);
    }  
 
  }
   
  boolean checkStruct(){
  
    boolean isOK=true;
	
	return isOK;
   
  }
  
  void initSeqObjAlignArray(List<SeqCompoAlignInfo> seqObjList, int seqTypeNum){  
	
    //int seqTypeNum=seqTypeInfo.seqTypeName.size();	
      ArrayList<Integer> seqAlignSStart;	
      ArrayList<Integer> seqAlignSEnd;  
    
      for(int s=0;s<seqObjList.size();s++){   
	        
			 seqAlignSStart=new ArrayList<Integer>();
			 seqAlignSEnd=new ArrayList<Integer>();
	         for(int i=0;i<seqTypeNum;i++){	           
				seqAlignSStart.add(-1);
				seqAlignSEnd.add(-1);
	         }
	         seqObjList.get(s).seqAlignSStart=seqAlignSStart;
	         seqObjList.get(s).seqAlignSEnd=seqAlignSEnd;
             seqAlignSStart=null;
			 seqAlignSEnd=null;		
		   		 
	  }	//for	
     
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
	int seqNum=0;
	try{    
		BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;
		String seqIdentifier;
		String seqLine;
		String seqName;
		String [] itemSplited;
		seqObjList=new ArrayList<SeqCompoAlignInfo>();
		seqNum=0;
		line = br.readLine();
	
		while(true){           
	       if (line == null) break;
	       line=line.trim();
		   if(line.indexOf(">")==0){
		     seqNum=seqNum+1;
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
			 perSeq.seqIdentifier=seqIdentifier;
			 perSeq.seqName=seqName;
			 perSeq.seqLength=seqLine.length();
			 perSeq.seq=seqLine;			 
			 seqObjList.add(perSeq);
			 perSeq=null;             
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
	int seqNum=0;
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
		seqNum=0;
	
		while(true){
           seqNum=seqNum+1;		   
           line = br.readLine();           	   
	       if (line == null) break;
		   line=line.trim();
           if(line.indexOf("@")!=0) {
		      System.out.println("Error in reading fastq");
			  break;
		   }			   
           seqIdentifier=line.substring(1,line.length());				 
		   itemSplited=seqIdentifier.split("\\s+");
		   seqName=itemSplited[0].trim();
		        
		   line=br.readLine();			 
		   if (line == null) break;
		   seqLine=line.trim();			 
		   seqLine=seqLine.replaceAll("N","n");
		   
		   line=br.readLine();				   
		   if (line == null) break;
		   line=line.trim();
		   if(line.indexOf("+")!=0) {
		      System.out.println("Error in reading FASTQ file");
			  break;
		   }
           
           line=br.readLine();			 
		   if (line == null) break;
		   seqQualityLine=line.trim();	
           
           perSeq=new SeqCompoAlignInfo();
           perSeq.seqIdentifier=seqIdentifier;
		   perSeq.seqName=seqName;
		   perSeq.seqLength=seqLine.length();
		   perSeq.seq=seqLine;				 
	       perSeq.seqQualityEncode=seqQualityLine;			 
		   seqObjList.add(perSeq);
		   perSeq=null;             		   
    
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
    SeqCompoAlignInfo perSeq=new SeqCompoAlignInfo();
	int seqNum=0;

	try{    

		BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;
		String seqIdentifier;
		String seqLine;
		String seqName;
		String [] itemSplited;
		seqObjList=new ArrayList<SeqCompoAlignInfo>();
		seqNum=0;
		line = br.readLine();				
		while(true){           
	       if (line == null) break;
	       line=line.trim();
		   if(line.indexOf(">")==0){
		     seqNum=seqNum+1;
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
			 seqLine=seqLine.replaceAll("N","n");
		  
			 perSeq=new SeqCompoAlignInfo();
			 perSeq.seqIdentifier=seqIdentifier;
			 perSeq.seqName=seqName;
			 perSeq.seqLength=seqLine.length();
			 perSeq.seq=seqLine;
			 seqObjList.add(perSeq);
			 perSeq=null;             
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
	int seqNum=0;
	
	try{    
       	BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;
		String seqIdentifier;	
		String seqLine;
		String seqName;		
		String seqQualityLine;
		String [] itemSplited;
	
		seqObjList=new ArrayList<SeqCompoAlignInfo>();
		seqNum=0;
			
		while(true){  
		   seqNum=seqNum+1;
		   
           line = br.readLine();			
	       if (line == null) break;
		   line=line.trim();
           if(line.indexOf("@")!=0) {
		      System.out.println("Error in reading fastq");
			  break;
		   }			   
           seqIdentifier=line.substring(1,line.length());				 
		   itemSplited=seqIdentifier.split("\\s+");
		   seqName=itemSplited[0].trim();
		      
		   line=br.readLine();			 
		   if (line == null) break;
		   seqLine=line.trim();			 
		   seqLine=seqLine.replaceAll("N","n");
		   
		   line=br.readLine();			  
		   if (line == null) break;
		   line=line.trim();
		   if(line.indexOf("+")!=0) {
		      System.out.println("Error in reading fastq");
			  break;
		   }	
           
           line=br.readLine();			 
		   if (line == null) break;
		   seqQualityLine=line.trim();
		   
		   perSeq=new SeqCompoAlignInfo();
		   perSeq.seqIdentifier=seqIdentifier;
		   perSeq.seqName=seqName;
		   perSeq.seqLength=seqLine.length();
		   perSeq.seq=seqLine;
		   perSeq.seqQualityEncode=seqQualityLine;
		   seqObjList.add(perSeq);
		   perSeq=null;	
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
 
 static boolean checkPairEndSeq(List<SeqCompoAlignInfo> seqObjList, 
		 List<SeqCompoAlignInfo> seqObjList2){
	
	 boolean isOK=true;
	 if(seqObjList.size()!=seqObjList2.size()){
		 System.out.println("Error: different seq num of pair-end seq");
		 return false;
	 }
	 //int nc=0;
	 for(int i=0;i<seqObjList.size();i++){
		if(!seqObjList.get(i).seqName.trim().equals(seqObjList2.get(i).seqName.trim())){ 
			isOK=false;
			System.out.println("Error: inconsistent seq name between pair-end seq");
			System.out.println("Forward: "+seqObjList.get(i).seqName
					+"\t Reverse:"+seqObjList2.get(i).seqName);
			break;
			//nc=nc+1;
			//System.out.println(seqObjList.get(i).seqName+"  "+seqObjList2.get(i).seqName);
		}
	 }
	 //System.out.println("inconsistent seq num:"+nc);
	 
	 return isOK;
  }
   
  
  void getBlastTarSeqOfRecognizer(List<SeqCompoAlignInfo> seqObjList, SeqRocket seqRocket, 
		  SeqCompoRecognizer recognizer, String outSeqFile){
  
	try{    
        int trimSEndIndex=-1;
		int trimLeftShift=0;
		int trimSEnd=0;
		for(int i=0;i<seqRocket.seqRecognizer.size();i++){
		  if(seqRocket.seqRecognizer.get(i).done 
				  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
		   trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
		   trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
		  }
		}
		
        int territoryLen=recognizer.territoryLen;
        boolean leftSubForBlast=recognizer.leftSubForBlast;    
		
		BufferedWriter writer=null;
        writer=new BufferedWriter(new FileWriter(outSeqFile));
		String nameLine;
		String seqLine="";
		String trimedSeqLine="";
	    int seqID=0;
		int seqLen=0;
    	if(trimSEndIndex>=0){
		  if(leftSubForBlast){
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
			  trimedSeqLine="nnnnnnnnnn";			
			  if(trimSEnd+territoryLen<seqLen)
				trimedSeqLine=seqLine.substring(trimSEnd,trimSEnd+territoryLen);
			  else if(trimSEnd<seqLen)
				trimedSeqLine=seqLine.substring(trimSEnd,seqLen);
				 
			  writer.write(trimedSeqLine);
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
			  trimedSeqLine="nnnnnnnnnn";	
			  if(trimSEnd<seqLen) trimedSeqLine=seqLine.substring(trimSEnd,seqLen);			 
			  
			  writer.write(trimedSeqLine);
			  writer.newLine();
			  writer.flush(); 
		    }		 
		  }
		}else{
		  if(leftSubForBlast){
		    for(int i=0;i<seqObjList.size();i++){	
			  seqID=i;	
			  nameLine=">"+seqID;
			  writer.write(nameLine);
			  writer.newLine();
			  //writer.flush(); 
			  seqLine=seqObjList.get(seqID).seq;
			  seqLen=seqLine.length();
				
			  if(territoryLen<seqLen)
				trimedSeqLine=seqLine.substring(0,territoryLen);
			  else 
				trimedSeqLine=seqLine.substring(0,seqLen);
				 
			  writer.write(trimedSeqLine);
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
  
  void getLeftSubSeq(List<SeqCompoAlignInfo> seqObj,int leftSubSeqLen,String outSeqFile){
	  		
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
  
  void recognizeExactBarcode(List<SeqCompoAlignInfo> seqObj, List<SeqRocket> rockets, 
		  String noMatchSeqFile, String inSeqFileFormat){
  		
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
		String barcodeSeq="";
		boolean isMatch=false;	
        int barIndex=0;
		String exactAlignedSeqFile="";
   		
		SeqCompoRecognizer barcodeInfo=null;
		ArrayList<BufferedWriter> writerList=new ArrayList<BufferedWriter> ();
		BufferedWriter writer0=null;
		for(int j=0;j<rockets.size();j++){
		  exactAlignedSeqFile=rockets.get(j).seqRecognizer.get(
			         rockets.get(j).seqTypeInfo.seqTypeName.indexOf(barcodeName)
			      ).exactAlignedSeqFile;
		  rockets.get(j).seqRecognizer.get(
				     rockets.get(j).seqTypeInfo.seqTypeName.indexOf(barcodeName)
				  ).exactAlignedSeqFile=exactAlignedSeqFile+"."+inSeqFileFormat;
		  writer0=new BufferedWriter(new FileWriter(exactAlignedSeqFile+"."+inSeqFileFormat));	
		  writerList.add(writer0);
          writer0=null;		  
		}
		BufferedWriter writer=null;
        writer=new BufferedWriter(new FileWriter(noMatchSeqFile));	
		if(saveAsFASTAFormat){
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;			
				isMatch=false;
				for(int j=0;j<rockets.size();j++){
				 barcodeInfo=rockets.get(j).seqRecognizer.get(
						 rockets.get(j).seqTypeInfo.seqTypeName.indexOf(barcodeName));
				 barcodeSeq=barcodeInfo.seq;
				 barIndex=seqLine.indexOf(barcodeSeq);			 
				 if(barIndex>=0 && barIndex<=barcodeInfo.maxExactStart-1){
					isMatch=true;
					writerList.get(j).write(">"+seqIdentifier);
					writerList.get(j).newLine();
					writerList.get(j).write(seqLine);
					writerList.get(j).newLine();
					writerList.get(j).flush(); 
				 } 
				}
				
				if(!isMatch){
				  writer.write(">"+seqIdentifier);
				  writer.newLine();
				  writer.write(seqLine);
				  writer.newLine();
				  writer.flush(); 
				}
			}
		}else if(saveAsFASTQFormat){
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;			
				seqQuality=seqObj.get(i).seqQualityEncode;
				
				isMatch=false;
				for(int j=0;j<rockets.size();j++){
				 barcodeInfo=rockets.get(j).seqRecognizer.get(
						 rockets.get(j).seqTypeInfo.seqTypeName.indexOf(barcodeName));
				 barcodeSeq=barcodeInfo.seq;
				 barIndex=seqLine.indexOf(barcodeSeq);			 
				 if(barIndex>=0 && barIndex<=barcodeInfo.maxExactStart-1){
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
			}
		}
        for(int j=0;j<writerList.size();j++){		
		   writerList.get(j).close();
		}
		writer.close();
		
		seqLine=null;
		seqIdentifier=null;		
		barcodeSeq=null;
		
	}catch(IOException e){
        System.out.println(e);
    }

 }

 void getRecognizedFASTAFile(List<SeqCompoAlignInfo> seqObj,SeqCompoRecognizer recognizer,
		 String seqType, SeqCompoType seqTypeInfo){ 
   	
    int alignSStartIndex=seqTypeInfo.seqTypeName.indexOf(seqType);	
	int alignSEndIndex=alignSStartIndex;	
	try{
		String seqLine;
		String seqIdentifier;			
		String trimedSeq="";
		
		boolean alignHTMLSave=recognizer.alignHTMLSave;
        boolean leftTrimSave=recognizer.leftTrimSave;
		boolean leftMaskSave=recognizer.leftMaskSave;
		boolean rightTrimSave=recognizer.rightTrimSave;
		boolean rightMaskSave=recognizer.rightMaskSave;
		int trimLeftShift=recognizer.trimLeftShift;		
		int trimRightShift=recognizer.trimRightShift;
        int trimIndex=0;		

		BufferedWriter writer0=null;	    
        writer0=new BufferedWriter(new FileWriter(recognizer.alignedFASTASeqFile));
		BufferedWriter writer1=null; 
		if(!leftTrimSave && !rightTrimSave){
		    for(int i=0;i<seqObj.size();i++){     
		     seqIdentifier=seqObj.get(i).seqIdentifier;
			 seqLine=seqObj.get(i).seq;				 
		
        	 writer0.write(">"+seqIdentifier);
			 writer0.newLine();
			 writer0.write(seqLine);
			 writer0.newLine();
			 writer0.flush(); 
			 
		    }		
		}else if(leftTrimSave && !rightTrimSave){
		    writer1=new BufferedWriter(new FileWriter(recognizer.alignTrimedSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;				
								
				writer0.write(">"+seqIdentifier);
				writer0.newLine();
				writer0.write(seqLine);
				writer0.newLine();
				writer0.flush(); 
				  
				trimIndex=seqObj.get(i).seqAlignSEnd.get(alignSEndIndex)-trimLeftShift;
				if(trimIndex<0) trimIndex=0;
				trimedSeq=seqLine.substring(trimIndex,seqLine.length());
				
				if (trimedSeq.length()==0) trimedSeq="nnnnnnnnnn";
				
				writer1.write(">"+seqIdentifier);
				writer1.newLine();
				writer1.write(trimedSeq);
				writer1.newLine();
				writer1.flush(); 	
                				
			}
			writer1.close();
		}else if(rightTrimSave && !leftTrimSave){
		    writer1=new BufferedWriter(new FileWriter(recognizer.alignTrimedSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;			
								
				writer0.write(">"+seqIdentifier);
				writer0.newLine();
				writer0.write(seqLine);
				writer0.newLine();
				writer0.flush(); 
				  
				trimIndex=seqObj.get(i).seqAlignSStart.get(alignSStartIndex)+trimRightShift;				
				if(seqLine.length()>=trimIndex)
				 trimedSeq=seqLine.substring(0,trimIndex);
				else if (seqLine.length()<trimIndex)
				 trimedSeq=seqLine.substring(0,seqLine.length());				 
				
				if (trimedSeq.length()==0) trimedSeq="nnnnnnnnnn";
				
				writer1.write(">"+seqIdentifier);
				writer1.newLine();
				writer1.write(trimedSeq);
				writer1.newLine();
				writer1.flush();
                				
			}
			writer1.close();
		} 
		
        writer0.close();
		
		if(alignHTMLSave){
           String line="";	
           String line1="";
           String line2="";
		   int sStart=0;
		   int sEnd=0;
			//int sStart0=0;
		   int sEnd0=0;	
            //int c=0;			
		   writer0=new BufferedWriter(new FileWriter(recognizer.alignedSeqHTMLFile));
		   writer0.write("<html><title></title><body>");
		   writer0.newLine();
		   writer0.flush(); 
		   for(int i=0;i<seqObj.size();i++){     
		     seqIdentifier=seqObj.get(i).seqIdentifier;
			 seqLine=seqObj.get(i).seq;					 
		
        	 writer0.write(">"+seqIdentifier+"<br>");
			 writer0.newLine();
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
                    		"<span style='color:"+seqTypeInfo.seqColor.get(k)+"'>";
				  }else{
				    sStart=sEnd0+1;
				    line1="<span style='color:"+seqTypeInfo.seqColor.get(k)+"'>";
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
			 if(sEnd0<seqLine.length())
			  line=line+seqLine.substring(sEnd0,seqLine.length());
			 line=line+"<br>";
			 
			 writer0.write(line);
			 writer0.newLine();
			 writer0.flush(); 
			 
		    }
            writer0.write("</body></html>");
			writer0.newLine();
			writer0.flush(); 
			
            writer0.close();			
		}
		
		String maskedSeq="";
		int maskIndex=0;
		StringBuilder str;
		String tmp="";
		if(leftMaskSave && !rightMaskSave){
		    writer1=new BufferedWriter(new FileWriter(recognizer.alignMaskedSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;						  
				maskIndex=seqObj.get(i).seqAlignSEnd.get(alignSEndIndex);
				if(maskIndex<0) maskIndex=0;
				
				str=new StringBuilder();	
				for (int x = 0; x < maskIndex; x++){
                  str.append("n");
                }
				tmp=str.toString();
				str=null;
                str=new StringBuilder(seqLine);				
				maskedSeq=str.replace(0,maskIndex,tmp).toString();
                str=null;		
				tmp=null;
				writer1.write(">"+seqIdentifier);
				writer1.newLine();
				writer1.write(maskedSeq);
				writer1.newLine();
				writer1.flush(); 	
                				
			}
			writer1.close();
		}else if(rightMaskSave && !leftMaskSave){
		    writer1=new BufferedWriter(new FileWriter(recognizer.alignMaskedSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
				seqIdentifier=seqObj.get(i).seqIdentifier;
				seqLine=seqObj.get(i).seq;				
			  
				maskIndex=seqObj.get(i).seqAlignSStart.get(alignSStartIndex);	
                if(maskIndex>seqLine.length()) maskIndex=seqLine.length();				
				
				str=new StringBuilder();	
				for (int x = 0; x < maskIndex; x++){
                  str.append("n");
                }
				tmp=str.toString();
				str=null;				
                str=new StringBuilder(seqLine);				
				maskedSeq=str.replace(maskIndex,seqLine.length(),tmp).toString();
                str=null;
				tmp=null;
				writer1.write(">"+seqIdentifier);
				writer1.newLine();
				writer1.write(maskedSeq);
				writer1.newLine();
				writer1.flush();
                				
			}
			writer1.close();
		} 
		
		
	}catch(IOException e){
        System.out.println(e);
    }	

 }
 
 void getRecognizedFASTQFile(List<SeqCompoAlignInfo> seqObj, String outFastqFile){    	

	try{
		String seqLine;
		String seqIdentifier;
        String seqQuality="";		
		
		BufferedWriter writer=null;	    
        writer=new BufferedWriter(new FileWriter(outFastqFile));
	
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
		
	}catch(IOException e){
        System.out.println(e);
    }	

 }
 
 List<SeqCompoAlignInfo> getRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, 
		 String seqType, SeqCompoType seqTypeInfo){
 
    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;		
	int alignSStartIndex=seqTypeInfo.seqTypeName.indexOf(seqType);		  
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
		 String seqType,SeqCompoType seqTypeInfo){
 
    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;	
	
	int alignSStartIndex=seqTypeInfo.seqTypeName.indexOf(seqType);		  
	//int alignSEndIndex=alignSStartIndex;
	
	for(int i=0; i<seqObjList.size(); i++){
	   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)==-1)
	     seqObj.add(seqObjList.get(i));
	
	}
	
	return seqObj;
 
 }
 
 List<SeqCompoAlignInfo> combineSeqObj(List<SeqCompoAlignInfo> seqObjList1, 
		 List<SeqCompoAlignInfo> seqObjList2){
 
    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo> ();
	
	for(int i=0; i<seqObjList1.size();i++){
	   seqObjList.add(seqObjList1.get(i));	   
	}
    
    for(int i=0; i<seqObjList2.size();i++){
	   seqObjList.add(seqObjList2.get(i));	   
	}
    
	return seqObjList;
 }
 
 List<Integer> getSeqIndex(List<SeqCompoAlignInfo> seqObj){
    ArrayList<Integer> seqIndex=new ArrayList<Integer>();
    for(int i=0;i<seqObj.size();i++){
	  seqIndex.add(seqObj.get(i).seqIndex);
	}
    return seqIndex;
 }
 
 List<SeqCompoAlignInfo> getSubSeqObj(List<SeqCompoAlignInfo> seqObj,
		 ArrayList<Integer> seqIndex){
 
    List<SeqCompoAlignInfo> subSeqObj=new ArrayList<SeqCompoAlignInfo>();
    for(int i=0;i<seqIndex.size();i++){
	  subSeqObj.add(seqObj.get(seqIndex.get(i)));
	}
    return subSeqObj;
 }
  
 List<SeqCompoAlignInfo> setSeqExactMatchPos(List<SeqCompoAlignInfo> seqObjList,
		 SeqCompoRecognizer recognizer, String seqType,SeqRocket seqRocket){
  
 	try{
		String seqLine;		
		int redSStartIndex=0;
        int redSStart=0;
        int redSEnd=0;		
        int redMaxStartIndex=recognizer.maxExactStart-1;
        if (redMaxStartIndex<0) redMaxStartIndex=0;
		
		int trimSEndIndex=-1;
		int trimLeftShift=0;
		int trimSEnd=0;		
		for(int i=0;i<seqRocket.seqRecognizer.size();i++){
		  if(seqRocket.seqRecognizer.get(i).done 
				  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
		    
			  trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
		      trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
		  
		  }
		}
		
		int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);		  
		int alignSEndIndex=alignSStartIndex;
		
		String recogniSeq=recognizer.seq;
		
        if(trimSEndIndex>=0){
		 for(int i=0;i<seqObjList.size();i++){

			seqLine=seqObjList.get(i).seq;				
			redSStartIndex=seqLine.indexOf(recogniSeq);
		    trimSEnd=seqObjList.get(i).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   	
			if(redSStartIndex>=trimSEnd && redSStartIndex<=trimSEnd+redMaxStartIndex){
			//if(redSStartIndex>=0){	 
				redSStart=redSStartIndex+1;
				redSEnd=redSStartIndex+recogniSeq.length();		
			}else{
			    redSStart=-1;
				redSEnd=-1;
			}
			
			seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,redSStart);
			seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,redSEnd);
			seqObjList.get(i).seqIndex=i;		
		 }
	   }else{
		 for(int i=0;i<seqObjList.size();i++){

			seqLine=seqObjList.get(i).seq;		
			
			redSStartIndex=seqLine.indexOf(recogniSeq);		   	
			if(redSStartIndex>=0 && redSStartIndex<=redMaxStartIndex){
			//if(redSStartIndex>=0){	 
				redSStart=redSStartIndex+1;
				redSEnd=redSStartIndex+recogniSeq.length();		
			}else{
			    redSStart=-1;
				redSEnd=-1;
			}
			
			seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,redSStart);
			seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,redSEnd);
			seqObjList.get(i).seqIndex=i;		
		 }
	   }		
    }catch(Exception e){
        System.out.println(e);
    }

	return seqObjList;

 }
   
 List<SeqCompoAlignInfo> getLeftSideBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
		 SeqCompoRecognizer recognizer, String seqType, SeqRocket seqRocket, 
		 String blastOutFile){
	
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
	    List<ArrayList <String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
		List <Integer> seqID=new ArrayList <Integer>();
		List <Integer> sStartList=new ArrayList <Integer>();
		List <Integer> sEndList=new ArrayList <Integer>();
		String readName="";
		int readID0=-1;
		int readID=0;
		
		for(int i=0; i<blastOut.size();i++){
		  readName=blastOut.get(i).get(blastInfoObj.colSName).trim();
		  readID=Integer.parseInt(readName);
		  if(readID>readID0){
				
			 alignLen=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colAlignLen).trim());
			 mismatchNum=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colMismatchNum).trim());
			 gapNum=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colGapNum).trim());
			 qStart=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colQStart).trim());
			 sStart=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSStart).trim());
			 sEnd=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSEnd).trim());	
	         /*
			 if(sStart>sEnd){
			   sEnd0=sStart;
			   sStart=sEnd;
			   sEnd=sEnd0;
			 }
			 */
	         maxMismatchNum=(int) Math.ceil(alignLen*maxMismatchRatio);
			 maxGapNum=(int) Math.ceil(alignLen*maxGapRatio);	 			
				
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
		  
		int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);	
		int alignSEndIndex=alignSStartIndex;	    
		
		int trimSEndIndex=-1;
		int trimLeftShift=0;
		int trimSEnd=0;
		for(int i=0;i<seqRocket.seqRecognizer.size();i++){
		  if(seqRocket.seqRecognizer.get(i).done 
				  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
			
			  trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
			  trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
		  
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
		
		seq=null;
		seqID=null;
		sStartList=null;
		sEndList=null;
		
		return seqObj;
	  
 } 
 
 
 List<ArrayList<SeqCompoAlignInfo>> getBaitBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
		 SeqCompoRecognizer recognizer, String seqType, SeqRocket seqRocket, 
		 String blastOutFile){
	
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
		
		//ArrayList<ArrayList <String>> blastOutList=new ArrayList<ArrayList <String>> ();
	    List<ArrayList <String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
		List <Integer> seqID=new ArrayList <Integer>();
		List <Integer> sStartList=new ArrayList <Integer>();
		List <Integer> sEndList=new ArrayList <Integer>();
		List <Integer> baitBrkQ=new ArrayList <Integer>();
		List <Integer> noSeqID=new ArrayList <Integer>();
		String readName="";
		int readID0=-1;
		int readID=0;
		
		for(int i=0; i<blastOut.size();i++){
		  readName=blastOut.get(i).get(blastInfoObj.colSName).trim();
		  readID=Integer.parseInt(readName);
		  if(readID>readID0){			
				alignLen=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colAlignLen).trim());
				mismatchNum=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colMismatchNum).trim());
				gapNum=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colGapNum).trim());
				qStart=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colQStart).trim());
				qEnd=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colQEnd).trim());
				sStart=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSStart).trim());
				sEnd=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSEnd).trim());	
	            if(sStart>sEnd){
				 sEnd0=sStart;
				 sStart=sEnd;
				 sEnd=sEnd0;
				}
	            maxMismatchNum=(int) Math.ceil(alignLen*maxMismatchRatio);
			    maxGapNum=(int) Math.ceil(alignLen*maxGapRatio);
				if(seqType.equals(baitName)){
					if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum 
							&& gapNum<=maxGapNum && qStart<=maxQStart && sStart<=maxSStart){
					  seqID.add(readID); 
					  sStartList.add(sStart);
					  sEndList.add(sEnd);				  
					  baitBrkQ.add(qEnd);
					  readID0=readID;
					 //blastOutList.add(blastOut.get(i));
					}
				}else if(seqType.equals(baitBrkName)){
				    baitBrkQPos=seqObjList.get(readID).baitBrkQPos;	
					maxSStart=recognizer.maxSStart-baitBrkQPos;
				    if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum && gapNum<=maxGapNum && qStart>baitBrkQPos && qStart<=maxQStart && sStart<=maxSStart){
					  seqID.add(readID); 
					  sStartList.add(sStart);
					  sEndList.add(sEnd);
					  readID0=readID;
					 //blastOutList.add(blastOut.get(i));
				    }
				
				}else if(seqType.equals(baitArmName)){
					if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum 
							&& gapNum<=maxGapNum && qStart<=maxQStart && sStart<=maxSStart){
					  seqID.add(readID); 
					  sStartList.add(sStart);
					  sEndList.add(sEnd);				  
					  readID0=readID;
					 //blastOutList.add(blastOut.get(i));
					}
				}
		  }
		}
	
		blastOut=null;		
		//blastOutList=null;
	  
		int alignSStartIndex=0;	
		alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);	
		int alignSEndIndex=alignSStartIndex;	    
		
		int trimSEndIndex=-1;
		int trimLeftShift=0;
		int trimSEnd=0;
		for(int i=0;i<seqRocket.seqRecognizer.size();i++){
		  if(seqRocket.seqRecognizer.get(i).done 
				  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
		   trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
		   trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
		  }
		}
		
		if(trimSEndIndex>=0){
			for(int i=0;i<seqID.size();i++){
			   seq =new SeqCompoAlignInfo();	
			   seq=seqObjList.get(seqID.get(i));
			   trimSEnd=seq.seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   		   
			   seq.seqAlignSStart.set(alignSStartIndex,trimSEnd+sStartList.get(i));
			   seq.seqAlignSEnd.set(alignSEndIndex,trimSEnd+sEndList.get(i));
			   if(seqType.equals(baitName)){
			    seq.baitBrkQPos=baitBrkQ.get(i);
			    seq.baitBrkSPos=trimSEnd+sEndList.get(i);
			   }
			   seqObj.add(seq);
			   seq=null;
			}
		}else{
			for(int i=0;i<seqID.size();i++){
			   seq =new SeqCompoAlignInfo();	
			   seq=seqObjList.get(seqID.get(i));
			   seq.seqAlignSStart.set(alignSStartIndex,sStartList.get(i));
			   seq.seqAlignSEnd.set(alignSEndIndex,sEndList.get(i));
			   if(seqType.equals(baitName)){
			    seq.baitBrkQPos=baitBrkQ.get(i);
			    seq.baitBrkSPos=sEndList.get(i);
			   }
			   seqObj.add(seq);
			   seq=null;
			}	
		}
		
		bait.add(seqObj);
		seqObj=null;
			
		for(int i=0; i<seqObjList.size();i++){
		   if(!seqID.contains(i)) {
		    noSeqID.add(i);
			//System.out.println(i);
		   }
		}
		seqID=null;
		sStartList=null;
		sEndList=null;
	    baitBrkQ=null;
		
		for(int i=0;i<noSeqID.size();i++){
		    seq =new SeqCompoAlignInfo();	
			seq=seqObjList.get(noSeqID.get(i));
			noSeqObj.add(seq);
			seq=null;
		}	
		
		bait.add(noSeqObj);
		noSeqObj=null;
		noSeqID=null;
		
		return bait;
	  
  }

   
  List<SeqCompoAlignInfo> shiftAlignPos(List<SeqCompoAlignInfo> seqObjList, 
		  String fromSeqType, String toSeqType, SeqCompoType seqTypeInfo){
	
    int fromSStartIndex=0;		  
	fromSStartIndex=seqTypeInfo.seqTypeName.indexOf(fromSeqType);
	int fromSEndIndex=fromSStartIndex;
    
    int toSStartIndex=0;		  
	toSStartIndex=seqTypeInfo.seqTypeName.indexOf(toSeqType);
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
	
	return seqObjList;
	  
 } 
  
 public List<SeqRocket> splitLaunchSingleEnd(String inSeqFile,int splitStep,String libSeqInfoFile,
		 String splitedSeqOut, String combinedSeqOut){
	 
	 List<String> splitedSeqFiles=new ArrayList<String>();
	//Check seq format, and then split seq into multiple subfiles................
	 System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
	 splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile, splitStep, splitedSeqOut);
	 
	 
	 List<SeqRocket> seqRockets=splitLaunchSingleEnd(splitedSeqFiles,libSeqInfoFile,combinedSeqOut);
	 
	 return seqRockets;
 }
 
 public List<SeqRocket> splitLaunchSingleEnd(List<String> splitedSeqFiles,String libSeqInfoFile,
		  String combinedSeqOut){

	 List<SeqRocket> seqRockets=new ArrayList<SeqRocket>();
	 
	 List<ArrayList<SeqRocket>> splitedRockets=new ArrayList<ArrayList<SeqRocket>>();
	 ArrayList<SeqRocket> rockets;
	 List<String> tmpFiles=new ArrayList<String>();
	 String outDir = null;
	 int s=1;
	 for(String file:splitedSeqFiles){
	    try {
		   outDir=file.substring(0,file.lastIndexOf("/"));
		   outDir=outDir+"/"+file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."));
		   FileOperate.newFolder(outDir);
		   
		   System.out.println("......Recognizing sequence for split "+s+" ......");

		   rockets=(ArrayList<SeqRocket>) launchSingleEnd(file,libSeqInfoFile,outDir);		  
	       if(rockets!=null && rockets.size()>0) splitedRockets.add(rockets);
		   
		   s++;
		   
		   tmpFiles.add(file);
		   for(SeqRocket rocket:rockets){
			  if(rocket.barcodeRecognizedSeqFile!=null)
			    tmpFiles.add(rocket.barcodeRecognizedSeqFile); 
			  if(rocket.recognizedSeqFile!=null)
			    tmpFiles.add(rocket.recognizedSeqFile); 
			  if(rocket.recognizedMaskSeqFile!=null)
			    tmpFiles.add(rocket.recognizedMaskSeqFile);
			  if(rocket.recognizedTrimSeqFile!=null)
			    tmpFiles.add(rocket.recognizedTrimSeqFile); 
		   }
		   rockets=null;
		}catch (Exception e) {
		// TODO Auto-generated catch block
		   e.printStackTrace();
		}
	 }
	   
	 if(splitedSeqFiles.size()>0){
		 String file=splitedSeqFiles.get(0);
		 if(combinedSeqOut==null) 
			 combinedSeqOut=file.substring(0,file.lastIndexOf("/"))+"/combined";
		 FileOperate.newFolder(combinedSeqOut);
		 System.out.println("......Combine splited sequences......");
		 List<String> fileList=new ArrayList<String>();
		 String outName="";
		 String seqFormat="fna";
	     for(int r=0;r<splitedRockets.get(0).size();r++){
	    	 SeqRocket rocket=new SeqRocket();
	    	 rocket=splitedRockets.get(0).get(r); 
			 
			 splitedSeqFiles=new ArrayList<String>();
			 for(int i=0;i<splitedRockets.size();i++){
				splitedSeqFiles.add(splitedRockets.get(i).get(r).recognizedSeqFile);
			 }
			 file=splitedSeqFiles.get(0);
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_combined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles,outName);
			 rocket.recognizedSeqFile=outName;
			 fileList.add(outName);
					 
			 splitedSeqFiles=new ArrayList<String>();			
			 for(int i=0;i<splitedRockets.size();i++){
				splitedSeqFiles.add(splitedRockets.get(i).get(r).recognizedMaskSeqFile);
			 }
			 file=splitedSeqFiles.get(0);
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_combined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles,outName);
			 rocket.recognizedMaskSeqFile=outName;
			 
	    	 splitedSeqFiles=new ArrayList<String>();	    	 
			 for(int i=0;i<splitedRockets.size();i++){
				splitedSeqFiles.add(splitedRockets.get(i).get(r).barcodeRecognizedSeqFile);
			 }
			 file=splitedSeqFiles.get(0);
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_combined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles,outName);
			 rocket.barcodeRecognizedSeqFile=outName;
			 
			 seqRockets.add(rocket);
			 rocket=null;
	     }
	     
	     String fileListOut=combinedSeqOut+"/fileList.txt";
		 FileOperate.saveList(fileList, null, fileListOut);

	 }
	 splitedRockets=null;
	 
	 setSeqRockets(seqRockets);	
	 
	 if(tmpFiles!=null){
		for(String tmpFile: tmpFiles){
		   FileOperate.delFile(tmpFile);
		}
	 }
	 
	 System.out.println("......Seq Recognization Done......");
	 
     return seqRockets;
 }
 
 public List<SeqRocketPair> splitLaunchPairEnd(String inSeqFile,String inSeqFile2,
		 String libSeqInfoFile,String libSeqInfoFile2,int splitStep,
		 String splitedSeqOut,String combinedSeqOut){
	 
	 List<String> splitedSeqFiles=new ArrayList<String>();
	 List<String> splitedSeqFiles2=new ArrayList<String>();
	 
	//Check seq format, and then split seq into multiple subfiles................	 
	 System.out.println("Total Seq Num(Pair-end Forward): "+SeqOperation.getSeqNum(inSeqFile));
	 splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile, splitStep, splitedSeqOut);
	 System.out.println("Total Seq Num(Pair-end Reverse): "+SeqOperation.getSeqNum(inSeqFile2));
	 splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2, splitStep, splitedSeqOut);

	 List<SeqRocketPair>seqPairRockets=splitLaunchPairEnd(splitedSeqFiles,splitedSeqFiles2,
			 libSeqInfoFile,libSeqInfoFile2,combinedSeqOut);	
	 
	 return seqPairRockets;

 }
 
 public List<SeqRocketPair> splitLaunchPairEnd(List<String>splitedSeqFiles,
		 List<String>splitedSeqFiles2,String libSeqInfoFile,String libSeqInfoFile2,
		 String combinedSeqOut){	 

	 List<SeqRocketPair> seqPairRockets=new ArrayList<SeqRocketPair>();

	 List<ArrayList<SeqRocketPair>> splitedRockets12=new ArrayList<ArrayList<SeqRocketPair>>();
	 ArrayList<SeqRocketPair> rockets12;
	 List<String> tmpFiles=new ArrayList<String>();
	 String outDir = null;
	 String file;
	 String file2;
	 for(int s=0;s<splitedSeqFiles.size();s++){
	    try {
	       file=splitedSeqFiles.get(s);
	       file2=splitedSeqFiles2.get(s);
	       outDir=file.substring(0,file.lastIndexOf("/"));
		   outDir=outDir+"/"+file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."));
		   FileOperate.newFolder(outDir);
		   
		   System.out.println("......Recognizing sequence for split "+(s+1)+" ......");
           
		   rockets12=(ArrayList<SeqRocketPair>) launchPairEnd(
				       file,file2,libSeqInfoFile,libSeqInfoFile2,outDir
				     );
		   if(rockets12!=null && rockets12.size()>0) splitedRockets12.add(rockets12);
		   
		   tmpFiles.add(file);
		   for(SeqRocketPair rocketPair:rockets12){
			  if(rocketPair.forward.barcodeRecognizedSeqFile!=null)
			    tmpFiles.add(rocketPair.forward.barcodeRecognizedSeqFile); 
			  if(rocketPair.forward.recognizedSeqFile!=null)
			    tmpFiles.add(rocketPair.forward.recognizedSeqFile); 
			  if(rocketPair.forward.recognizedMaskSeqFile!=null)
			    tmpFiles.add(rocketPair.forward.recognizedMaskSeqFile);
			  if(rocketPair.forward.recognizedTrimSeqFile!=null)
			    tmpFiles.add(rocketPair.forward.recognizedTrimSeqFile); 
			  
			  if(rocketPair.reverse.barcodeRecognizedSeqFile!=null)
				tmpFiles.add(rocketPair.reverse.barcodeRecognizedSeqFile); 
			  if(rocketPair.reverse.recognizedSeqFile!=null)
				tmpFiles.add(rocketPair.reverse.recognizedSeqFile); 
			  if(rocketPair.reverse.recognizedMaskSeqFile!=null)
			    tmpFiles.add(rocketPair.reverse.recognizedMaskSeqFile);
			  if(rocketPair.reverse.recognizedTrimSeqFile!=null)
			    tmpFiles.add(rocketPair.reverse.recognizedTrimSeqFile); 
		   }
		   rockets12=null;
	
		}catch (Exception e) {
		// TODO Auto-generated catch block
		   e.printStackTrace();
		}
	 }
	   
	 if(splitedRockets12.size()>0){
		 file=splitedSeqFiles.get(0);
		 if(combinedSeqOut==null) 
			 combinedSeqOut=file.substring(0,file.lastIndexOf("/"))+"/combined";
		 FileOperate.newFolder(combinedSeqOut);
		 System.out.println("......Combine splited sequences......");
		 List<String> fileList_forward=new ArrayList<String>();	
		 List<String> fileList_reverse=new ArrayList<String>();		
		 String outName="";
		 String seqFormat="fna";
		 for(int r=0;r<splitedRockets12.get(0).size();r++){
			 SeqRocketPair rocketPair=new SeqRocketPair();
			 rocketPair.forward=splitedRockets12.get(0).get(r).forward;
			 rocketPair.reverse=splitedRockets12.get(0).get(r).reverse;
			 
			 // for recognized Seq
			 List<String> splitedSeqFiles_F=new ArrayList<String>();
			 List<String> splitedSeqFiles_R=new ArrayList<String>();
			 for(int i=0;i<splitedRockets12.size();i++){
				splitedSeqFiles_F.add(
				   splitedRockets12.get(i).get(r).forward.recognizedSeqFile
				);
				splitedSeqFiles_R.add(
				   splitedRockets12.get(i).get(r).reverse.recognizedSeqFile
				);
			 }
			 file=splitedSeqFiles_F.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ForwardCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
			 rocketPair.forward.recognizedSeqFile=outName;
			 fileList_forward.add(outName);
			 
			 file=splitedSeqFiles_R.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ReverseCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName);
			 rocketPair.reverse.recognizedSeqFile=outName;
			 fileList_reverse.add(outName);
			 
			 // for recognized Seq masked
			 splitedSeqFiles_F=new ArrayList<String>();
			 splitedSeqFiles_R=new ArrayList<String>();
			 for(int i=0;i<splitedRockets12.size();i++){
				splitedSeqFiles_F.add(
				   splitedRockets12.get(i).get(r).forward.recognizedMaskSeqFile
				);
				splitedSeqFiles_R.add(
				   splitedRockets12.get(i).get(r).reverse.recognizedMaskSeqFile
				);
			 }
			 file=splitedSeqFiles_F.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ForwardCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
			 rocketPair.forward.recognizedMaskSeqFile=outName;
			 
			 file=splitedSeqFiles_R.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ReverseCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName);
			 rocketPair.reverse.recognizedMaskSeqFile=outName;
			 
			 //for barcode-recognized Seq
			 splitedSeqFiles_F=new ArrayList<String>();
			 splitedSeqFiles_R=new ArrayList<String>();
			 for(int i=0;i<splitedRockets12.size();i++){
				splitedSeqFiles_F.add(
				   splitedRockets12.get(i).get(r).forward.barcodeRecognizedSeqFile
				);
				splitedSeqFiles_R.add(
				   splitedRockets12.get(i).get(r).reverse.barcodeRecognizedSeqFile
				);
			 }
			 file=splitedSeqFiles_F.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ForwardCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_F,outName);
			 rocketPair.forward.barcodeRecognizedSeqFile=outName;
			 
			 file=splitedSeqFiles_R.get(0);	
			 seqFormat=FileOperate.getFileFormat(file);
			 outName=combinedSeqOut+"/"
			         +file.substring(file.lastIndexOf("/")+1,file.lastIndexOf("."))
		             +"_ReverseCombined."+seqFormat;
			 SeqOperation.combineSeqFile(splitedSeqFiles_R,outName);
			 rocketPair.reverse.barcodeRecognizedSeqFile=outName;
			 
			 seqPairRockets.add(rocketPair);
			 rocketPair=null;
			 splitedSeqFiles_F=null;
			 splitedSeqFiles_R=null;
	     }	
		 
		 String fileListOut_forward=combinedSeqOut+"/fileList_forward.txt";
		 FileOperate.saveList(fileList_forward, null, fileListOut_forward);
		 String fileListOut_reverse=combinedSeqOut+"/fileList_reverse.txt";
		 FileOperate.saveList(fileList_reverse, null, fileListOut_reverse);
		 fileList_forward=null;
		 fileList_reverse=null;

	 }
	 splitedRockets12=null;
	 
	 setSeqPairRockets(seqPairRockets);
	 
	 if(tmpFiles!=null){
		for(String tmpFile: tmpFiles){
		   FileOperate.delFile(tmpFile);
		}
	 }
	 
	 System.out.println("......Seq Recognization Done......");
	 
	 return seqPairRockets;
	 
 }

 public List<SeqRocket> launchSingleEnd(String inSeqFile,String libSeqInfoFile,
		 String seqOutDir) {	    
		
		FileOperate.newFolder(tmpDir+"/recognizer");
		FileOperate.newFolder(tmpDir+"/recognizer/barcode");
		FileOperate.newFolder(tmpDir+"/recognizer/primer");
		FileOperate.newFolder(tmpDir+"/recognizer/bait");
		
		tmpFiles=new  ArrayList<String>();
				
		if(seqOutDir==null){
		   seqOutDir=inSeqFile.substring(0, inSeqFile.lastIndexOf("/"))+"/RecognizedSeq";	
		}
		FileOperate.newFolder(seqOutDir);
		seqOutDir=seqOutDir+"/";

		// forward seq(Single-end)
		List<SeqRocket> rockets =buildRockets(libSeqInfoFile);		
				
		List<SeqCompoAlignInfo> seqObjList=null; 
		//List<SeqCompoAlignInfo> payLoadSeqObj;
				
		//Setting and Checkpoint for pair-end forward Seq................
		System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
		seqObjList=checkSeq(inSeqFile);
		System.out.println("Checked Seq Num: "+seqObjList.size());
		String inSeqFileFormat=FileOperate.getFileFormat(inSeqFile);
		
		//Recognizing barcode ...............................................	
		
		String barNoExactSeqFile=tmpDir+"/AllBarNoExactMatchSeq_forward."+inSeqFileFormat;
        String barNoExactSeqLeftSubFile=barNoExactSeqFile.substring(
        		                          0,barNoExactSeqFile.lastIndexOf(".")
        	                            )+".leftsub.fna";       
		
        tmpFiles.add(barNoExactSeqFile);
        tmpFiles.add(barNoExactSeqLeftSubFile);
        
		List<SeqCompoAlignInfo> barNoExactSeqObj=getBarNoExactSeqObj(seqObjList,rockets,
				inSeqFileFormat, barNoExactSeqFile,barNoExactSeqLeftSubFile);		
		seqObjList=null;
		
		for(int i=0;i<rockets.size();i++){		  
		  
			launchSeqRocket(rockets.get(i), barNoExactSeqObj,
					barNoExactSeqLeftSubFile,inSeqFileFormat, seqOutDir);			
			

		}// for Rocket
		barNoExactSeqObj=null;    
	    
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
	      FileOperate.delFile(tmpFile);
		}
	    
	    setSeqRockets(rockets);
	    
	    return rockets;
	   
 }
 
 public List<SeqRocketPair> launchPairEnd(String inSeqFile,String inSeqFile2,
		 String libSeqInfoFile,String libSeqInfoFile2,String seqOutDir) {	    
		
		FileOperate.newFolder(tmpDir+"/recognizer");
		FileOperate.newFolder(tmpDir+"/recognizer/barcode");
		FileOperate.newFolder(tmpDir+"/recognizer/primer");
		FileOperate.newFolder(tmpDir+"/recognizer/bait");
		
		tmpFiles=new  ArrayList<String>();
				
		if(seqOutDir==null){
		   seqOutDir=inSeqFile.substring(0, inSeqFile.lastIndexOf("/"))+"/RecognizedSeq";	
		}
		FileOperate.newFolder(seqOutDir);
		seqOutDir=seqOutDir+"/";
		
		
		List<SeqCompoAlignInfo> seqObjList=null; //pair-end forward
		List<SeqCompoAlignInfo> seqObjList2=null; //pair-end reversed
		//List<SeqCompoAlignInfo> payLoadSeqObj;
				
		//Setting and Checkpoint for pair-end forward Seq................
		System.out.println("Total Seq Num(pair-end forward): "
		                    +SeqOperation.getSeqNum(inSeqFile)
		                  );
		String inSeqFileFormat=FileOperate.getFileFormat(inSeqFile);
		seqObjList=checkSeq(inSeqFile);
		System.out.println("Checked Seq Num(pair-end forward): "+seqObjList.size());
		
		// Setting and Checkpoint for pair-end reversed Seq................
		System.out.println("Total Seq Num(pair-end reverse): "+SeqOperation.getSeqNum(inSeqFile2));
		String inSeqFileFormat2=FileOperate.getFileFormat(inSeqFile2);
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
		
		List<SeqRocket> rockets =buildRockets(libSeqInfoFile);	// for forward seq(Pair-end)
		//Recognizing barcode for pair-end forward Seq...............................................	
		String barNoExactSeqFile=tmpDir+"/AllBarNoExactMatchSeq_forward."+inSeqFileFormat;
        String barNoExactSeqLeftSubFile=barNoExactSeqFile.substring(
        		                          0,barNoExactSeqFile.lastIndexOf(".")
        	                            )+".leftsub.fna";  
        tmpFiles.add(barNoExactSeqFile);
        tmpFiles.add(barNoExactSeqLeftSubFile);
			
		List<SeqCompoAlignInfo> barNoExactSeqObj=getBarNoExactSeqObj(seqObjList,rockets,
				inSeqFileFormat, barNoExactSeqFile,barNoExactSeqLeftSubFile);
		seqObjList=null;		

		for(int i=0;i<rockets.size();i++){		  
		  
			launchSeqRocket(rockets.get(i), barNoExactSeqObj,
					barNoExactSeqLeftSubFile,inSeqFileFormat,seqOutDir);
			
			String forwardSeqFile=rockets.get(i).recognizedSeqFile;
			String forwardSeqFormat=FileOperate.getFileFormat(forwardSeqFile);
			String forwardMaskSeqFile=rockets.get(i).recognizedMaskSeqFile;
			String forwardMaskSeqFormat=FileOperate.getFileFormat(forwardMaskSeqFile);
			String forwardSeqNameFile=tmpDir+"/"+forwardSeqFile.substring(
					 forwardSeqFile.lastIndexOf("/")+1,forwardSeqFile.lastIndexOf(".")
			    	)+".seqName";
		    tmpFiles.add(forwardSeqNameFile);

			SeqOperation.extratSeqName(forwardSeqFile,forwardSeqNameFile);
			String reverseSeqFile=forwardSeqFile.substring(
					                0,forwardSeqFile.lastIndexOf(".")
					              )+".Reverse."+inSeqFileFormat2;
			SeqOperation.getSubSeq(inSeqFile2,forwardSeqNameFile,0,reverseSeqFile);
			seqObjList2=checkSeq(reverseSeqFile);
			
			//Recognizing barcode for pair-end reversed Seq...............................................	
			String barNoExactSeqFile2=tmpDir+"/AllBarNoExactMatchSeq_reverse."+inSeqFileFormat2;
	        String barNoExactSeqLeftSubFile2=barNoExactSeqFile2.substring(
	        		   0,barNoExactSeqFile2.lastIndexOf(".")
	        	   )+".leftsub.fna";
	        tmpFiles.add(barNoExactSeqFile2);
	        tmpFiles.add(barNoExactSeqLeftSubFile2);
						
	        List<SeqRocket> rockets2 =buildRockets(libSeqInfoFile2); // for reversed seq(Pair-end)
			List<SeqCompoAlignInfo> barNoExactSeqObj2=getBarNoExactSeqObj(seqObjList2,
				rockets2,inSeqFileFormat2, barNoExactSeqFile2,barNoExactSeqLeftSubFile2);
			seqObjList2=null;

			for(int j=0;j<rockets2.size();j++){	
			   try {  
				    
				    String rocket12Tag=rockets.get(i).rocketName+"-"
				                       +rockets2.get(j).rocketName;
				    
				    // get reverse-recognized result
				    SeqRocket rocket_reverse=rockets2.get(j);
				    rocket_reverse.rocketName=rocket12Tag+"_Reverse";					
					launchSeqRocket(rocket_reverse, barNoExactSeqObj2,
							barNoExactSeqLeftSubFile2, inSeqFileFormat2, seqOutDir);
							
					// get forward seq based on reverse-recognized result
					SeqRocket rocket_forward=(SeqRocket) rockets.get(i).clone();
					rocket_forward.rocketName=rocket12Tag+"_Forward";
					String reverseSeqOut=rocket_reverse.recognizedSeqFile;							
					String reverseSeqNameFile=tmpDir+"/"+reverseSeqOut.substring(
						reverseSeqOut.lastIndexOf("/")+1,reverseSeqOut.lastIndexOf(".")
				      )+".seqName";	
			        tmpFiles.add(reverseSeqNameFile);			       
					SeqOperation.extratSeqName(reverseSeqOut,reverseSeqNameFile);
					
					String forwardSeqOut=reverseSeqOut.substring(
							          0,reverseSeqOut.lastIndexOf(".")
							       )+".Forward."+forwardSeqFormat;				
					SeqOperation.getSubSeq(forwardSeqFile,reverseSeqNameFile,0,forwardSeqOut);
					rocket_forward.recognizedSeqFile=forwardSeqOut;
                    
					String reverseMaskSeqOut=rocket_reverse.recognizedMaskSeqFile;	
					String forwardMaskSeqOut=reverseMaskSeqOut.substring(
							         0,reverseMaskSeqOut.lastIndexOf(".")
							       )+".Forward."+forwardMaskSeqFormat;
					SeqOperation.getSubSeq(forwardMaskSeqFile,reverseSeqNameFile,0,forwardMaskSeqOut);
					rocket_forward.recognizedMaskSeqFile=forwardMaskSeqOut;			
				
					SeqRocketPair rocketPair=new SeqRocketPair();
					rocketPair.forward=rocket_forward;
					rocketPair.reverse=rocket_reverse;
					rockets12.add(rocketPair);
					
					rocketPair=null;
					rocket12Tag=null;
					reverseSeqOut=null;
					reverseSeqNameFile=null;
					forwardSeqOut=null;
					forwardMaskSeqOut=null;
				}catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}			
			}
			rockets2=null;
			barNoExactSeqObj2=null;
			barNoExactSeqFile2=null;
			barNoExactSeqLeftSubFile2=null;
			reverseSeqFile=null;
			forwardSeqFile=null;
			forwardSeqFormat=null;
			forwardMaskSeqFile=null;
			forwardMaskSeqFormat=null;
			forwardSeqNameFile=null;
					
		}// for Rocket
		barNoExactSeqObj=null;    
	    
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
	      FileOperate.delFile(tmpFile);
		}
 	    
	    setSeqPairRockets(rockets12);	    
	 
	    //sSystem.out.println("Done!");
	    return rockets12;
	  
 }
 
 public List<SeqCompo> getLibExpSeqCompo(String libSeqInfoFile){
	    
	    List<SeqCompo> expSeqCompoList = new ArrayList<SeqCompo>();
	    
	    String expName=null;
	    String barName=null;
		String barSeq=null;	
		String primerSeq=null;	
		String primerContSeq=null;	
		String baitSeq=null;
		String baitArmSeq=null;
		String freqCutterSeq=null;
		
		boolean isExpNameOK=false;	
		boolean isPrimerOK=false;
	    int expIdx=0;
	    List<String> expNameList=new ArrayList<String>();
		try{    
		   BufferedReader br;              
		   br = new BufferedReader(new FileReader(libSeqInfoFile));
		   String line;
		   String [] itemSplited;
		   String attrName;
		   String attrValue;
		   line = br.readLine();		
		   while(true){ 			   
			  if(line == null) break;
			  if(line.trim().indexOf(">")==0){			     
				 line=br.readLine();			
				 if (line == null) break;
				 expName=null;
				 barName=null;
				 barSeq=null;	
				 primerSeq=null;	
				 primerContSeq=null;	
				 baitSeq=null;
				 baitArmSeq=null;
				 isExpNameOK=false;
				 isPrimerOK=false;
				 expIdx=expIdx+1;
				 while(line.indexOf(">")<0 && line.indexOf("=")>0){							 
					itemSplited=line.split("=");
					if(itemSplited.length>=2 && itemSplited[0]!=null 
							&& itemSplited[1]!=null){
					  
						attrName=itemSplited[0].trim();
						attrValue=itemSplited[1].trim();
					
						if(attrName.equalsIgnoreCase("ExperimentName")){
						   expName=attrValue.trim();
						   if(expName==null){
							  System.out.println("!!! Warning: You have empty 'ExperimentName'");
							  System.out.println("We automatically use it's index as experiment name.");
							  expName="Exp_"+expIdx; 
						   }
						   
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
							 barName=attrValue.trim();							
						}else if(attrName.equalsIgnoreCase("BarcodeSeq")){
							 barSeq=attrValue.trim();					
						}else if(attrName.equalsIgnoreCase("PrimerSeq")){
							 primerSeq=attrValue.trim();
							 isPrimerOK=true;
						}else if(attrName.equalsIgnoreCase("PrimerContSeq")){
							 primerContSeq=attrValue.trim();	
						}else if(attrName.equalsIgnoreCase("BaitTerritorySeq")){
							 baitSeq=attrValue.trim();	
						}else if(attrName.equalsIgnoreCase("BaitTerritoryArmSeq")){
							 baitArmSeq=attrValue.trim();	
						}else if(attrName.equalsIgnoreCase("FrequentCutterSeq")){
							 freqCutterSeq=attrValue.trim();	
						}						
				
						line=br.readLine(); 
						    
					}else{
						line=br.readLine();						
					}
					 
					if (line == null) break;
					    
				 }
				 
				 if(isPrimerOK){
					 if(!isExpNameOK){
						System.out.println("!!! Warning: You didn't correctly configure 'ExperimentName' for No."+expIdx);
						System.out.println("We automatically use it's index as experiment name.");
						expName="Exp_"+expIdx; 
						expNameList.add(expName);
					 }
					 SeqCompo seqComp = new SeqCompo();
					 seqComp.isActive=true;
					 seqComp.name=expName;
					 seqComp.barName=barName;
					 seqComp.barSeq=barSeq;
					 seqComp.primerSeq=primerSeq;
					 seqComp.primerContSeq=primerContSeq;
					 seqComp.baitTerritorySeq=baitSeq;
					 seqComp.baitTerritoryArmSeq=baitArmSeq;
					 seqComp.freqCutterSeq=freqCutterSeq;
				
					 expSeqCompoList.add(seqComp);
					 seqComp=null;							
				 }else{
					 System.out.println("!!! Warning: You didn't correctly configure 'PrimerSeq' for No."+expIdx);
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
  
 List<SeqRocket> buildRockets(String libSeqInfoFile){
	 
		List<SeqCompo> libExpSeqCompList=getLibExpSeqCompo(libSeqInfoFile);
			
		SeqRocket libSeqRocket;		
		SeqCompoRecognizer barcode;
	    SeqCompoRecognizer primer;
	    SeqCompoRecognizer primerCont;
		SeqCompoRecognizer bait;
		SeqCompoRecognizer baitBrk;
		SeqCompoRecognizer baitArm;
		
		String barSeq=null;	
		String primerSeq=null;	
		String primerContSeq=null;	
		String baitSeq=null;
		String baitArmSeq=null;		
		String rocketName="";	
		
		List<SeqRocket> rockets =new ArrayList<SeqRocket>();
		for(int i=0;i<libExpSeqCompList.size();i++){
		    barSeq=libExpSeqCompList.get(i).barSeq;	
		    primerSeq=libExpSeqCompList.get(i).primerSeq;	
		    primerContSeq=libExpSeqCompList.get(i).primerContSeq;	
		    baitSeq=libExpSeqCompList.get(i).baitTerritorySeq;
		    baitArmSeq=libExpSeqCompList.get(i).baitTerritoryArmSeq;	
		    
			rocketName=libExpSeqCompList.get(i).name;	
			if(libExpSeqCompList.get(i).barName!=null)
			  rocketName=rocketName+"."+libExpSeqCompList.get(i).barName;
			
			libSeqRocket=new SeqRocket();
			libSeqRocket.rocketName=rocketName;
			libSeqRocket.isActive=libExpSeqCompList.get(i).isActive;
		
			libSeqRocket.seqRecognizer=new ArrayList<SeqCompoRecognizer> ();
			libSeqRocket.seqTypeInfo=new SeqCompoType();
			libSeqRocket.seqTypeInfo.seqTypeName=new ArrayList<String> ();	
			libSeqRocket.seqTypeInfo.seqColor=new ArrayList<String> ();
			libSeqRocket.seqComponent=libExpSeqCompList.get(i);
			if(primerSeq!=null && primerSeq.length()>0){
			  if(primerSeq.length()<minRedPrimerLen){
			    System.out.println("The length of your primer is less than "+minRedPrimerLen+"!");
			    minRedPrimerLen=primerSeq.length();
			  }
			  barcode=new SeqCompoRecognizer();
			  primer=new SeqCompoRecognizer();
			  libSeqRocket.seqTypeInfo.seqTypeName.add(barcodeName);
			  libSeqRocket.seqTypeInfo.seqTypeName.add(primerName);
			  libSeqRocket.seqTypeInfo.seqColor.add(barcodeColor);
			  libSeqRocket.seqTypeInfo.seqColor.add(primerColor);
			  barcode.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(barcodeName);	
			  primer.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(primerName);
			  barcode.rawSeq=barSeq;	
			  primer.rawSeq=primerSeq;
			  barcode.seqName=rocketName+".Barcode";	
			  primer.seqName=rocketName+".Primer";	
			  libSeqRocket.seqRecognizer.add(barcode);
			  libSeqRocket.seqRecognizer.add(primer);  
			}
			
			if(baitSeq!=null && baitSeq.length()>0){
			  primerCont=new SeqCompoRecognizer();		
			  libSeqRocket.seqTypeInfo.seqTypeName.add(primerContName);
			  libSeqRocket.seqTypeInfo.seqColor.add(primerContColor);
			  primerCont.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(primerContName);
			  primerCont.rawSeq=primerContSeq;	
			  primerCont.seqName=rocketName+".PrimerCont";
	          libSeqRocket.seqRecognizer.add(primerCont);
	          bait=new SeqCompoRecognizer();		  
	          libSeqRocket.seqTypeInfo.seqTypeName.add(baitName);		  
	          libSeqRocket.seqTypeInfo.seqColor.add(baitColor);		  
			  bait.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(baitName);		  	
			  bait.rawSeq=baitSeq;		  
			  bait.seqName=rocketName+".Bait";		  
			  libSeqRocket.seqRecognizer.add(bait);
			  baitBrk=new SeqCompoRecognizer();		  
			  libSeqRocket.seqTypeInfo.seqTypeName.add(baitBrkName);		  
			  libSeqRocket.seqTypeInfo.seqColor.add(baitBrkColor);		  
			  baitBrk.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(baitBrkName);		  	
			  baitBrk.rawSeq=baitSeq;		  
			  baitBrk.seqName=rocketName+".BaitBrk";		  
			  libSeqRocket.seqRecognizer.add(baitBrk);
			  if(baitArmSeq!=null && baitArmSeq.length()>0){
				  baitArm=new SeqCompoRecognizer();		  
				  libSeqRocket.seqTypeInfo.seqTypeName.add(baitArmName);		  
				  libSeqRocket.seqTypeInfo.seqColor.add(baitArmColor);		  
				  baitArm.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(baitArmName);		  	
				  baitArm.rawSeq=baitArmSeq;		  
				  baitArm.seqName=rocketName+".BaitArm";		  
				  libSeqRocket.seqRecognizer.add(baitArm);			  
			  }

			}else if(primerContSeq!=null && primerContSeq.length()>0){
			  primerCont=new SeqCompoRecognizer();
			  libSeqRocket.seqTypeInfo.seqTypeName.add(primerContName);
			  libSeqRocket.seqTypeInfo.seqColor.add(primerContColor);
			  primerCont.index=libSeqRocket.seqTypeInfo.seqTypeName.indexOf(primerContName);	
			  primerCont.rawSeq=primerContSeq;
			  primerCont.seqName=rocketName+".primerCont";		  
			  libSeqRocket.seqRecognizer.add(primerCont);		     		  
			}			
			
			rockets.add(libSeqRocket);
			libSeqRocket=null;
			barcode=null;
			primer=null;
			primerCont=null;
		    bait=null;
			baitBrk=null;
			baitArm=null;
		
		}
		
		libExpSeqCompList=null;
		
		configRecognizer(rockets);
		
		return rockets;
 }
  
 List<SeqCompoAlignInfo> getBarNoExactSeqObj(List<SeqCompoAlignInfo>seqObjList,
		List<SeqRocket> rockets,String inSeqFileFormat,
		String barNoExactSeqFile,String leftSubBarNoExactSeqFile){
	  
		//Recognizing barcode ...............................................		
		
		recognizeExactBarcode(seqObjList,rockets,barNoExactSeqFile,inSeqFileFormat);	
		
		for(int f=0;f<rockets.size();f++){
			tmpFiles.add(
			  rockets.get(f).seqRecognizer.get(
				    rockets.get(f).seqTypeInfo.seqTypeName.indexOf(barcodeName)
			  ).exactAlignedSeqFile
			);
		}		
					
		List<SeqCompoAlignInfo> barNoExactSeqObj=getSeqObj(barNoExactSeqFile);
	
		getLeftSubSeq(barNoExactSeqObj,barTerritoryLen,leftSubBarNoExactSeqFile);	
		
		tmpFiles.add(barNoExactSeqFile);
		tmpFiles.add(leftSubBarNoExactSeqFile);
		
		return barNoExactSeqObj;
 }
  
 void launchSeqRocket(SeqRocket seqRocket, List<SeqCompoAlignInfo> barNoExactSeqObj, 
		String barNoExactSeqLeftSubFile, String inSeqFileFormat,String seqOutDir){		
		
		  String barBlastCMD="";
		  String tarBlastCMD="";
		  String baitBlastCMD="";
		  String blastCMD="";	
		 
		  String tarSeqFile="";
		  String tarBaitSeqFile="";
		  String querySeqFile="";
		  String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		  String blastOutFile=tmpDir+"/SeqRocket.BlastOut."+timeStamp+".txt";

		  int seqNum=0;
		  int blastWordSize=7;
		  String blastTask="blastn-short";		
		  
		  List<SeqCompoAlignInfo> tarSeqObj = null;
		  List<SeqCompoAlignInfo> tarExactSeqObj;
		  List<SeqCompoAlignInfo> tarNoExactSeqObj;
		  List<SeqCompoAlignInfo> tarBlastSeqObj;	
		  List<ArrayList<SeqCompoAlignInfo>> baitObj;
		  List<Integer> blastWordSizeList;
		  String seqType="";
		  String lastSeqType="";
	      boolean doAlignment=false;		
		  boolean doBarBlast=true;
		  boolean leftTrim=false;
		  
		  seqNum=SeqOperation.getSeqNum(barNoExactSeqLeftSubFile);		
		  barBlastCMD="-evalue 10000 -max_target_seqs "+seqNum;	
		  if(seqNum==0) doBarBlast=false;
		  String doName="";
		  String saveName="";
		  SeqCompoRecognizer recognizer;
		  SeqCompoRecognizer primerCont;
		  
		  if(seqRocket.isActive){
			tarSeqObj=new ArrayList<SeqCompoAlignInfo> ();
			doName=seqRocket.rocketName+".";
			doAlignment=true;
			if(doAlignment){
			    initSeqObjAlignArray(barNoExactSeqObj,
			    		seqRocket.seqTypeInfo.seqTypeName.size());
			}			
			for(int k=0;k<seqRocket.seqRecognizer.size();k++){
			  
			   if(!seqRocket.seqTypeInfo.seqTypeName.contains(barcodeName)){
				  System.out.println("Warning: "
			           +"You did not set barcode or primer sequence for experiment '"
				       +doName+"', so no sequence extraction for this experiment!!!");
					  
				  break;
			   }				  
			   
			   recognizer=seqRocket.seqRecognizer.get(k);	
			   seqType=seqRocket.seqTypeInfo.seqTypeName.get(recognizer.index);	  
			   System.out.println(recognizer.seqName+" Recognizing......");
			   doName=doName+seqType+"-";	
			  
			   if(seqType.equals(barcodeName)){             		  
				  
				  if((inSeqFileFormat.equalsIgnoreCase("fastq") 
						  || inSeqFileFormat.equalsIgnoreCase("fnq") 
						  || inSeqFileFormat.equalsIgnoreCase("fq"))){
				     
					  recognizer.saveSeqAsFASTQFormat=true;
					  recognizer.saveSeqAsFASTAFormat=false;
				  }				  
				 			  
				  tarExactSeqObj=getSeqObj(recognizer.exactAlignedSeqFile);
	              
				  //doAlignment=true;
				  if(doAlignment){  
				      initSeqObjAlignArray(tarExactSeqObj,
						  seqRocket.seqTypeInfo.seqTypeName.size());
				  }
				  
				  tarExactSeqObj=setSeqExactMatchPos(tarExactSeqObj,recognizer,
	            		  seqType,seqRocket); 
				  tarBlastSeqObj=new ArrayList<SeqCompoAlignInfo> (); 			  
				  if(doBarBlast){
					  querySeqFile=recognizer.tagSeqFile;
					  blastWordSize=recognizer.blastWordSize;
					  blastTask=recognizer.blastTask;
					  blastCMD="blastn -task "+blastTask+" "
					          +barBlastCMD+" -word_size="+blastWordSize
					          +" -query "+querySeqFile
					          +" -subject "+barNoExactSeqLeftSubFile
					          +" -out "+blastOutFile+" -outfmt 6";
					  //System.out.println(blastCMD+"......");
					  SeqOperation.runBLAST2Seq(blastCMD);			  
					  tarBlastSeqObj=getLeftSideBLASTSeq(barNoExactSeqObj,recognizer,
							  seqType,seqRocket,blastOutFile);		  
					  FileOperate.delFile(blastOutFile);				 		  
				  }  
	              lastSeqType=seqType;				  
	              tarSeqObj=combineSeqObj(tarExactSeqObj,tarBlastSeqObj);
	              System.out.println(seqRocket.rocketName+">>>Got "
				      +seqType+" Recognized Seq: "+tarSeqObj.size()
				      +" (Exact:"+tarExactSeqObj.size()
				      +", NoExact:"+tarBlastSeqObj.size()+")");			  
				 
			   }else if(seqType.equals(baitName) || seqType.equals(baitBrkName) 
					  || seqType.equals(baitArmName)){
				  
				  if((inSeqFileFormat.equalsIgnoreCase("fastq") 
							  || inSeqFileFormat.equalsIgnoreCase("fnq") 
							  || inSeqFileFormat.equalsIgnoreCase("fq"))){
					     
					  recognizer.saveSeqAsFASTQFormat=true;
						
				  }	
				  tarExactSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				  tarNoExactSeqObj=tarSeqObj;
				  tarSeqObj=null;
				  tarSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				  primerCont=null;	
				  
	              blastWordSizeList=new ArrayList<Integer> ();
				  if(seqType.equals(baitArmName) || seqType.equals(baitBrkName)){
				     if(recognizer.seqLength>=2000) blastWordSizeList.add(128);
					 if(recognizer.seqLength>=1000) blastWordSizeList.add(64);
					 if(recognizer.seqLength>=500) blastWordSizeList.add(48);
					 if(recognizer.seqLength>=50) blastWordSizeList.add(21);
				     blastWordSizeList.add(11); 
				  }else if(seqType.equals(baitName)){
					  if(recognizer.seqLength>=500) blastWordSizeList.add(48);
					  if(recognizer.seqLength>=100) blastWordSizeList.add(32);
					  blastWordSizeList.add(21);
					  blastWordSizeList.add(16);
					  blastWordSizeList.add(11);
					  blastWordSizeList.add(7);				  
					  blastWordSizeList.add(4);
					  primerCont=seqRocket.seqRecognizer.get(
							  seqRocket.seqTypeInfo.seqTypeName.indexOf(primerContName));	
				  } 			  
				  for(int c=0;c<blastWordSizeList.size();c++){
					  recognizer.blastWordSize=blastWordSizeList.get(c);
					  recognizer.minAlignLen=recognizer.blastWordSize;
					  if(seqType.equals(baitName) && blastWordSizeList.get(c)<=primerCont.blastWordSize){
						 recognizer.leftSubForBlast=true;
						 recognizer.territoryLen=primerCont.territoryLen;
						 recognizer.minAlignLen=primerCont.minAlignLen;
						 recognizer.blastWordSize=primerCont.blastWordSize;
						 recognizer.tagSeqFile=primerCont.tagSeqFile;
					  }
					  if(recognizer.leftSubForBlast && leftTrim){
						tarBaitSeqFile=tmpDir+"/"+seqType+".LeftTrim.LeftSub.ForBLAST."+c+".fna";	   
					  }else if(recognizer.leftSubForBlast){
						tarBaitSeqFile=tmpDir+"/"+seqType+".LeftSub.ForBLAST."+c+".fna";	   
					  }if(leftTrim){
					    tarBaitSeqFile=tmpDir+"/"+seqType+".LeftTrim.ForBLAST."+c+".fna";	   
					  }else{
						tarBaitSeqFile=tmpDir+"/"+seqType+".ForBLAST."+c+".fna";	   
					  }
					  getBlastTarSeqOfRecognizer(tarNoExactSeqObj,seqRocket,
							  recognizer,tarBaitSeqFile);
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
					  baitBlastCMD="blastn -task "+blastTask
							  +" -word_size="+blastWordSize+" -evalue 10000";
					  seqNum=SeqOperation.getSeqNum(tarBaitSeqFile);		
					  baitBlastCMD=baitBlastCMD+" -max_target_seqs "+seqNum;         	  
					  querySeqFile=recognizer.tagSeqFile;
					  blastCMD=baitBlastCMD+" -query "+querySeqFile+" -subject "
					          +tarBaitSeqFile+" -out "+blastOutFile+" -outfmt 6";	
					  tarBlastSeqObj=new ArrayList<SeqCompoAlignInfo> ();
					  baitObj=new ArrayList<ArrayList<SeqCompoAlignInfo>>();
					  if(seqNum>0){
						//System.out.println(blastCMD+"......");
						SeqOperation.runBLAST2Seq(blastCMD);		
						baitObj=getBaitBLASTSeq(tarNoExactSeqObj,recognizer,
								seqType,seqRocket,blastOutFile);
						FileOperate.delFile(blastOutFile);
						tarBlastSeqObj=baitObj.get(0);   //seqObj having BLAST alignment,meaning those seq are in bait breaksite or bait arm region.
						tarNoExactSeqObj=baitObj.get(1); //seqObj not having BLAST alignment, meaning those seq are not in bait breaksite or bait arm region.
						if(tarBlastSeqObj.size()>0) 
							tarSeqObj=combineSeqObj(tarSeqObj,tarBlastSeqObj);						
					  }else{
						//System.out.println(" No Sequence for "+recognizer.seqName+" Blast!");
						break;
					  }	  
					  baitObj=null;
					  	
				  }// for blastWordSize	
				  if(tarNoExactSeqObj.size()>0){			    
					  tarSeqObj=combineSeqObj(tarSeqObj,tarNoExactSeqObj);
					  //Set alignArray for seqObj without BLAST alignment in order to show that this seq junction goes to out of bait region.
					  tarSeqObj=shiftAlignPos(tarSeqObj,lastSeqType,seqType,seqRocket.seqTypeInfo);
				  }
				  System.out.println(seqRocket.rocketName+">>>Seq harbouring "+seqType
					   +" portion/total recognized seq: "
					   +(tarSeqObj.size()-tarNoExactSeqObj.size()) +"/"+tarSeqObj.size());
				  tarNoExactSeqObj=null;
	              lastSeqType=seqType;			  
			  
			   }else{
				  if((inSeqFileFormat.equalsIgnoreCase("fastq") 
							  || inSeqFileFormat.equalsIgnoreCase("fnq") 
							  || inSeqFileFormat.equalsIgnoreCase("fq"))){
					     
					  recognizer.saveSeqAsFASTQFormat=true;
						
				  }	
					  
				  tarSeqObj=setSeqExactMatchPos(tarSeqObj,recognizer,seqType,seqRocket);  
				  tarExactSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				  tarNoExactSeqObj=new ArrayList<SeqCompoAlignInfo> ();	  
				  tarExactSeqObj=getRecognizedSeq(tarSeqObj,seqType,seqRocket.seqTypeInfo);
				  tarNoExactSeqObj=getNoRecognizedSeq(tarSeqObj,seqType,seqRocket.seqTypeInfo);
				  tarSeqObj=null;			 
				  if(recognizer.leftSubForBlast && leftTrim){
					tarSeqFile=tmpDir+"/"+seqType+".Trim.leftsub.ForBLAST.fna";	   
				  }else if(recognizer.leftSubForBlast){
					tarSeqFile=tmpDir+"/"+seqType+".LeftSub.ForBLAST.fna";	   
				  }else if(leftTrim){
				    tarSeqFile=tmpDir+"/"+seqType+".Trim.ForBLAST.fna";	
				  }else{
					tarSeqFile=tmpDir+"/"+seqType+".ForBLAST.fna";	    
				  }	
				  getBlastTarSeqOfRecognizer(tarNoExactSeqObj,seqRocket,recognizer,tarSeqFile);	
	              tmpFiles.add(tarSeqFile);		  
				  blastWordSize=recognizer.blastWordSize;
	              blastTask=recognizer.blastTask;			  
				  tarBlastCMD="blastn -task "+blastTask
						  +" -word_size="+blastWordSize+" -evalue 10000";
				  seqNum=SeqOperation.getSeqNum(tarSeqFile);		
				  tarBlastCMD=tarBlastCMD+" -max_target_seqs "+seqNum;         	  
				  querySeqFile=recognizer.tagSeqFile;
				  blastCMD=tarBlastCMD+" -query "+querySeqFile
						  +" -subject "+tarSeqFile+" -out "+blastOutFile+" -outfmt 6";
				  tarBlastSeqObj=new ArrayList<SeqCompoAlignInfo> ();
				  if(seqNum>0){
					//System.out.println(blastCMD+"......");
				    SeqOperation.runBLAST2Seq(blastCMD);		
				    tarBlastSeqObj=getLeftSideBLASTSeq(tarNoExactSeqObj,recognizer,
							seqType,seqRocket,blastOutFile);          
					FileOperate.delFile(blastOutFile);			
				  }else{
					//System.out.println(" No Sequence for "+recognizer.seqName+" Blast!");
				  }	 
				  tarNoExactSeqObj=null;	 
				  tarSeqObj=combineSeqObj(tarExactSeqObj,tarBlastSeqObj);
	              lastSeqType=seqType;
				  System.out.println(seqRocket.rocketName+">>> Got "+seqType
					  +" Recognized Seq: "+tarSeqObj.size()
					  +" (Exact:"+tarExactSeqObj.size()+", NoExact:"+tarBlastSeqObj.size()+")");
				  
			   }
			  
			   tarExactSeqObj=null;
			   tarBlastSeqObj=null;
			  
			   if(recognizer.saveRecognizedSeq){
				   saveName=doName.substring(0,doName.length()-1);
				   //saveName=seqRocket.rocketName;
				   //saveName=recognizer.seqName;
				   if(seqRocket.seqRecognizer.get(k).saveSeqAsFASTAFormat){
				      seqRocket.seqRecognizer.get(k).alignedFASTASeqFile
	                              =seqOutDir+saveName+".AlignedSeq.fna"; 
				      seqRocket.seqRecognizer.get(k).alignMaskedSeqFile
	                              =seqOutDir+saveName+".AlignedSeqMasked.fna";
				      seqRocket.seqRecognizer.get(k).alignTrimedSeqFile
				                  =seqOutDir+saveName+".AlignedSeqTrimed.fna";				   
				      seqRocket.seqRecognizer.get(k).alignedSeqHTMLFile
	                              =seqOutDir+saveName+".AlignedSeq.html";				   
				      getRecognizedFASTAFile(tarSeqObj,seqRocket.seqRecognizer.get(k),
						   seqType,seqRocket.seqTypeInfo);			
				   
				      seqRocket.recognizedMaskSeqFile
				                  =seqRocket.seqRecognizer.get(k).alignMaskedSeqFile;				   
				      seqRocket.recognizedSeqFile
	                              =seqRocket.seqRecognizer.get(k).alignedFASTASeqFile;
				      if(seqType.equals(barcodeName)){
					       seqRocket.barcodeRecognizedSeqFile
                                  =seqRocket.seqRecognizer.get(k).alignedFASTASeqFile;
				      }
				   }
				   
				   if(seqRocket.seqRecognizer.get(k).saveSeqAsFASTQFormat){
					  seqRocket.seqRecognizer.get(k).alignedFASTQSeqFile
					              =seqOutDir+saveName+".AlignedSeq.fastq"; 	
					  
					  getRecognizedFASTQFile(tarSeqObj,
							   seqRocket.seqRecognizer.get(k).alignedFASTQSeqFile);	
					  
					  if(seqRocket.recognizedSeqFile!=null) 
						  FileOperate.delFile(seqRocket.recognizedSeqFile);
					  
					  seqRocket.recognizedSeqFile
                              =seqRocket.seqRecognizer.get(k).alignedFASTQSeqFile;
					  
					  if(seqType.equals(barcodeName)){
						 seqRocket.barcodeRecognizedSeqFile
	                          =seqRocket.seqRecognizer.get(k).alignedFASTQSeqFile;
					  }
				   }
				   
				   System.out.println(seqRocket.rocketName +">>>Saved Recognized "+saveName+" Seq!");	  
			   }
			  
			   seqRocket.seqRecognizer.get(k).done=true;		
			  
			   if(recognizer.leftShiftForNext){
				 leftTrim=true;	
			   }
			   recognizer=null;	  
			} // for seqRecognizer
			
			tarSeqObj=null; //release memory space

		  }// is rocket active?         
  }
 
}
  
  
