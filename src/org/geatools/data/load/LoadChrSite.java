package org.geatools.data.load;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.collections.Predicate;
import org.geatools.NumberCheck;
import org.geatools.data.structure.ChrInfo;
import org.geatools.data.structure.ChrSite;
import org.geatools.operation.FileOperation;

public class LoadChrSite {
	 static int dataStartLine=0;	 
	 static int chrColIdx=0;
	 static int chrStartColIdx=1;
	 static int chrEndColIdx=2;	
	 static int nameColIdx=3;
	 static int scoreColIdx=4;
	 static int chrStrandColIdx=5;
	 
	 static int genTSSUpStreamLen=5000;
	 static int genEndDownStreamLen=1000;
	 static String chrRegex=LoadChrInfo.getChrRegex();
	 static String chr;
    

	 public static List <ChrSite> getChrSites(String chrSiteFile){
		   
		    if(chrSiteFile==null || !new File(chrSiteFile).exists()) return null;
		    
		    List <ChrSite> chrSiteList=new ArrayList <ChrSite> ();
		    ChrSite chrSite;
			List <ArrayList <String>> chrSites=FileOperation.getMatrixFromFile(chrSiteFile);

		    String chrom="";
	        String strand="";
	        String score;
			int index=-1;
			for(int i=0;i<chrSites.size();i++){
			   chrom=chrSites.get(i).get(0).trim();
			   if(chrom.matches(chrRegex)){
				   index=index+1;				   
				   chrSite= new ChrSite ();	
				   chrSite.index=index;
				   chrSite.chr=chrSites.get(i).get(chrColIdx).trim();
				   chrSite.chrStart=Integer.parseInt(chrSites.get(i).get(chrStartColIdx).trim());
				   chrSite.chrEnd=Integer.parseInt(chrSites.get(i).get(chrEndColIdx).trim());
				   strand="";
				   if(chrSites.get(i).size()<(chrStrandColIdx+1) 
						   || chrSites.get(i).get(chrStrandColIdx)==null){
					  chrSite.strand=".";
				   }else {					   
					  strand=chrSites.get(i).get(chrStrandColIdx).trim();
					  if(strand.equals("+") || strand.equals("-"))
					     chrSite.strand=strand;
					  else
					     chrSite.strand=".";  						   
				   }
				   
				   if(chrSites.get(i).size()>nameColIdx){
				     chrSite.name=chrSites.get(i).get(nameColIdx).trim();
				   }else{
					 chrSite.name=".";
				   }
				   
				   if(chrSites.get(i).size()>scoreColIdx){
					  score=chrSites.get(i).get(scoreColIdx).trim();
					  if(NumberCheck.isNumeric(score)) 
						  chrSite.score=Double.parseDouble(score);
					  else 
						  chrSite.score=null;						  
				   }else{
					  chrSite.score=null;
				   }
				   				  	
			       chrSiteList.add(chrSite);
				   
			       chrSite=null;
			   }
			}
			
			chrSites=null;			
		
			return chrSiteList;
			
	 }
	 
	 public static List<ChrSite> getChrSites(List<ChrInfo> chrInfoList){
		 
		 List<ChrSite> chrSiteList=new ArrayList<ChrSite>();		
		 int index=0;
		 for(int i=0;i<chrInfoList.size();i++){	
		    //Collections.sort(chrInfoList.get(i).siteList,new ChrSite.CompSiteStart(false));	
			for(ChrSite site:chrInfoList.get(i).siteList){	
			  site.index=index;
			  chrSiteList.add(site);
			  index=index+1;
			}
		 }
		 
		 return chrSiteList;
	 }

	 
	 public static void setChrSites(String chrSiteFile,List<ChrInfo> chrInfoList){		   

			if(chrSiteFile==null || !new File(chrSiteFile).exists()) return;
			
			ChrSite chrSite;
			List <ArrayList <String>> chrSites=FileOperation.getMatrixFromFile(chrSiteFile);
			
			String chrom="";
			String strand="";
			String score;
			for(ChrInfo chrInfo:chrInfoList){
			  chrInfo.siteList=new ArrayList<ChrSite>();
			}
			
			int index=-1;
			for(int i=0;i<chrSites.size();i++){
				chrom=chrSites.get(i).get(0).trim();
				if(chrom.matches(chrRegex)){
				  index=index+1;					   
				  chrSite= new ChrSite();	
				  chrSite.index=index;
				  chrSite.chr=chrSites.get(i).get(chrColIdx).trim();
				  chrSite.chrStart=Integer.parseInt(chrSites.get(i).get(chrStartColIdx).trim());
				  chrSite.chrEnd=Integer.parseInt(chrSites.get(i).get(chrEndColIdx).trim());
				  strand="";
				  if(chrSites.get(i).size()<(chrStrandColIdx+1) 
						   || chrSites.get(i).get(chrStrandColIdx)==null){
					  chrSite.strand="+";
				  }else {					   
					  strand=chrSites.get(i).get(chrStrandColIdx).trim();
					  if(strand.equals("+") || strand.equals("-")) chrSite.strand=strand;
					  else chrSite.strand="+";  						   
				  }
				  
				  if(chrSites.get(i).size()>nameColIdx){
				     chrSite.name=chrSites.get(i).get(nameColIdx).trim();
				  }else{
					 chrSite.name=".";
				  }
				  
				  if(chrSites.get(i).size()>scoreColIdx){
					  score=chrSites.get(i).get(scoreColIdx).trim();
					  if(NumberCheck.isNumeric(score)) 
						  chrSite.score=Double.parseDouble(score);
					  else 
						  chrSite.score=null;						  
				  }else{
					  chrSite.score=null;
				  }
				  
				  chr=chrSite.chr;
				  ChrInfo chrInfo=(ChrInfo) CollectionUtils.find(chrInfoList,new Predicate(){
				        public boolean evaluate(Object o) {
				           return chr.equals(((ChrInfo) o).name);
				        }
				  });
				  
				  if(chrInfo!=null) chrInfo.siteList.add(chrSite);
				  
				  chrSite=null;				
		
				}
			}			
			
			for(int i=0;i<chrInfoList.size();i++){	
			   Collections.sort(chrInfoList.get(i).siteList,new ChrSite.CompSiteStart(false));	
			}
    }	
	 
	public static void saveChrSites(List<ChrSite> chrSiteList,String outFile){
	     
		 if(chrSiteList==null) return;
		 
		 try{	         
			 String text="";
	         BufferedWriter writer=null;
	         writer=new BufferedWriter(new FileWriter(outFile));

	         for(ChrSite chrSite:chrSiteList){
			    text=chrSite.chr+"\t"+chrSite.chrStart+"\t"+chrSite.chrEnd+"\t"
	                 +chrSite.name+"\t"+chrSite.score+"\t"+chrSite.strand;
	             
	            writer.write(text);
	            writer.newLine();
	            writer.flush();    				  
	         } 
				
			 writer.close();
			 writer=null;
	         text=null;
	     }catch(IOException e){
	          System.out.println(e);
	     }
	}
	
	public static void saveChrSiteSeqAsFASTA(List<ChrSite> chrSiteList,String outFile){
	     
		 if(chrSiteList==null) return;
		 
		 try{	         
			 String name="";		
	         BufferedWriter writer=null;
	         writer=new BufferedWriter(new FileWriter(outFile));

	         for(ChrSite chrSite:chrSiteList){
			    name=">"+chrSite.name+" "+chrSite.chr+":"+chrSite.chrStart+"-"+chrSite.chrEnd+" "
	                 +"strand("+chrSite.strand+")";	             
	            writer.write(name);
	            writer.newLine();
	    		writer.write(chrSite.seq);
	    		writer.newLine();
	            writer.flush();    				  
	         } 				
			 writer.close();
			 writer=null;
	        
	     }catch(IOException e){
	          System.out.println(e);
	     }
	 }

}
