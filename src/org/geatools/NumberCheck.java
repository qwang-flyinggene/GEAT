package org.geatools;

public class NumberCheck {
	 
		  /*
	 public static boolean isNumeric(String str){
	     return str.matches("^(?:(?:\\-{1})?\\d+(?:\\.{1}\\d+)?)$");
	 }
         */
	 public static boolean isNumeric(String str){  
	      
		 if(str==null || str.trim()==null) return false;
		 try{  
	    	  Double.parseDouble(str);     	        
	     }catch(NumberFormatException nfe){  
	    	  return false;  
	     }   
	   
	     return true;  
	 }
	 
	 public static boolean isPositiveNumeric(String str){  
		 
		 if(str==null || str.trim()==null) return false;
		 
		 boolean isPositive=false;
	     try{  
	    	  double d = Double.parseDouble(str); 
	    	  if(d>=0 ) isPositive=true;
	     }catch(NumberFormatException nfe){  
	    	  return false;  
	     }   
	   
	     return isPositive;  
	 }
	 
	 public static boolean isNegativeNumeric(String str){  
		  if(str==null || str.trim()==null) return false;
		  boolean isNegative=false;
	      try{  
	    	  double d = Double.parseDouble(str); 
	    	  if(d<=0 ) isNegative=true;
	      }catch(NumberFormatException nfe){  
	    	  return false;  
	      }   
	   
	      return isNegative;  
	 }
	 
	 public static boolean isPercentile(String str){
		  if(str==null || str.trim()==null) return false;
	      boolean isPercentile=false;
	      try{  
	        double d = Double.parseDouble(str);  
	        if(d>=0 && d<=1) isPercentile=true;	        
	      }catch(NumberFormatException nfe){  
	        return false;  
	      } 
	      
	      return isPercentile;  
	 }
	 
	 public static boolean isPValue(String str){
		  if(str==null || str.trim()==null) return false;
	      boolean isPValue=false;
	      try{  
	        double d = Double.parseDouble(str);  
	        if(d>=0 && d<=1) isPValue=true;	        
	      }catch(NumberFormatException nfe){  
	        return false;  
	      } 
	      
	      return isPValue;  
	  }
	 
	  public static boolean isInteger(String strNum) {
		  boolean ret = true;
		  try {
		     Integer.parseInt(strNum);
		  }catch (NumberFormatException e) {
		     ret = false;
		  }
		 return ret;
	  }
	    
	  public static boolean isPositiveInteger(String str){
		  if(str==null || str.trim()==null) return false;
	        boolean isPositive=false;
	        try{  
	          double d = Integer.parseInt(str);  
	          if(d>=0) isPositive=true;
	          
	        }catch(NumberFormatException nfe){  
	          return false;  
	        } 
	        
	        return isPositive;  
	  }
	  
	  public static boolean isNegativeInteger(String str){
		    if(str==null || str.trim()==null) return false;
	        boolean isNegative=false;
	        try{  
	          double d = Integer.parseInt(str);  
	          if(d<=0) isNegative=true;
	          
	        }catch(NumberFormatException nfe){  
	          return false;  
	        } 
	        
	        return isNegative;  
	  }
	  
	  public static double getMaxValue(double[][] data){
		  
		  double maxValue=Double.MIN_VALUE;
		  for(int i=0;i<data.length;i++){
			  for(int j=0;j<data[i].length;j++){
				  if(maxValue<data[i][j]) maxValue=data[i][j];
			  }			  
		  }
		  
		  return maxValue;
		  
	  }
	  
     public static double getMinValue(double[][] data){
		  
		  double minValue=Double.MAX_VALUE;
		  for(int i=0;i<data.length;i++){
			  for(int j=0;j<data[i].length;j++){
				  if(minValue>data[i][j]) minValue=data[i][j];
			  }			  
		  }
		  
		  return minValue;
		  
	  }
     
     public static double getDigitLengthOfNum(double number){
   	  double digitNum=0;
   	  
   	  if(number>0){
   	    digitNum= Math.floor(Math.log10(number) + 1);
   	  }else if(number<0){
   		digitNum= Math.floor(Math.log10(-number) + 1); 
   	  }
   	  
   	  return digitNum;
     }
}