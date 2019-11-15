;; This file contains a number of functions and programs written by
;; OTHER programmers throughout the years, with either small
;; modifications by JZG or put here because significantly different
;; versions are floating around with which TAP is not compatible.

;; Written by John Johnson:
;; hline
;; vline
;; sharpcorners
;; sigfig

;; Written by JAV
;; which

;; Written by  W. Landsman, modified by many:
;; reacol (tap_readcol)

;; Han Wen
;; strreplace

;; D. Lindler
;; strnumber (with W. Landsman, MRG RITSS)
;; remchar

;; H. Freudenreich, STX, Landsman
;; robust_sigma


pro remchar,st,char	;Remove character
;+
; NAME:
;	REMCHAR
; PURPOSE:
;	Remove all appearances of character (char) from string (st)
;
; CALLING SEQUENCE:
;	REMCHAR, ST, CHAR
;
; INPUT-OUTPUT:
;	ST  - String from which character will be removed, scalar or vector  
; INPUT:
;	CHAR- Single character to be removed from string or all elements of a
;		string array 
;
; EXAMPLE:
;	If a = 'a,b,c,d,e,f,g' then 
;
;	IDL> remchar,a, ','
;
;      will give a = 'abcdefg'
;
; REVISIONS HISTORY
;	Written D. Lindler October 1986
;	Test if empty string needs to be returned   W. Landsman  Feb 1991
;	Work on string arrays    W. Landsman   August 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
                                 ;Convert string to byte
 if N_params() LT 2 then begin
     print,'Syntax - REMCHAR, string, character'
     return
 endif

 bchar = byte(char) & bchar = bchar[0]          ;Convert character to byte

 for i = 0,N_elements(st)-1 do  begin

 bst = byte(st[i])
 good = where( bst NE bchar, Ngood)
 if Ngood GT 0 then st[i] = string(bst[good]) else st[i] = ''

 endfor
 return
 end


FUNCTION  ROBUST_SIGMA,Y, ZERO=REF
;
;+
; NAME:
;ROBUST_SIGMA 
;
; PURPOSE:
;Calculate a resistant estimate of the dispersion of a distribution.
; EXPLANATION:
;For an uncontaminated distribution, this is identical to the standard
;deviation.
;
; CALLING SEQUENCE:
;result = ROBUST_SIGMA( Y, [ /ZERO ] )
;
; INPUT:
;Y = Vector of quantity for which the dispersion is to be calculated
;
; OPTIONAL INPUT KEYWORD:
;/ZERO - if set, the dispersion is calculated w.r.t. 0.0 rather than the
;central value of the vector. If Y is a vector of residuals, this
;should be set.
;
; OUTPUT:
;ROBUST_SIGMA returns the dispersion. In case of failure, returns
;value of -1.0
;
; PROCEDURE:
;Use the median absolute deviation as the initial estimate, then weight
;points using Tukey's Biweight. See, for example, "Understanding Robust
;and Exploratory Data Analysis," by Hoaglin, Mosteller and Tukey, John
;Wiley & Sons, 1983.
;
; REVSION HISTORY:
;H. Freudenreich, STX, 8/90
;       Replace MED() call with MEDIAN(/EVEN)  W. Landsman   December 2001
;
;-  
EPS = 1.0E-20
 IF KEYWORD_SET(REF) THEN Y0=0. ELSE Y0  = MEDIAN(Y,/EVEN) 
; First, the median absolute deviation MAD about the median:  
MAD = MEDIAN( ABS(Y-Y0), /EVEN )/0.6745 
; If the MAD=0, try the MEAN absolute deviation:

 IF MAD LT EPS THEN MAD = AVG( ABS(Y-Y0) )/.80
 IF MAD LT EPS THEN RETURN, 0.0 ; Now the biweighted value:
 U   = (Y-Y0)/(6.*MAD)
 UU  = U*U
 Q   = WHERE(UU LE 1.0, COUNT)
 IF COUNT LT 3 THEN BEGIN
;   PRINT,'ROBUST_SIGMA: This distribution is TOO WEIRD! Returning -1'
   SIGGMA = -1.
   RETURN,SIGGMA
 ENDIF  
 
NUMERATOR = TOTAL( (Y[Q]-Y0)^2 * (1-UU[Q])^4 )
 N     = N_ELEMENTS(Y)
 DEN1  = TOTAL( (1.-UU[Q])*(1.-5.*UU[Q]) )
 SIGGMA = N*NUMERATOR/(DEN1*(DEN1-1.))  
IF SIGGMA GT 0. THEN RETURN, SQRT(SIGGMA) ELSE RETURN, 0.  
END

;+
; NAME:
;        SIGFIG
;
;
; PURPOSE:
;        Accept a scalar numerical value or an array of numbers and
;        return the numbers as strings with the specified number of
;        significant figures.
;
; CALLING SEQUENCE:
;        RESULT = SigFig(Number, Nfig [, /SCIENTIFIC, /PLUSSES, /NUMERICAL)
;
; INPUTS:
;        Number - Scalar or array of numerical values (float, double, int)
;        Nfig   - Number of signficant figures desired in the output
;
; OUTPUTS:
;        String representation of the input with the specified number
;        of signficant figures.
;
; KEYWORD PARAMTERS:
;        /SCIENTIFIC - return the numbers in scientific notation
;        /PLUSSES    - Include plus signs for positive numbers 
;        /NUMERICAL  - Return numerical, rather than string, values
;
; EXAMPLE:
;        IDL> print, sigfig(-0.0001234, 2)      
;        -0.00012
;        IDL> print, sigfig(1.234, 1)
;        1.
;        IDL> print, sigfig(1234, 1) 
;        1000
;        IDL> print, sigfig(-0.0001234, 2, /sci)
;        -1.2e-4
;        IDL> print, sigfig(1234, 2, /plus)
;        +1200
;        IDL> print, sigfig(1234, 2, /plus, /sci)
;        +1.2e+3
;
; MODIFICATION HISTORY:
; Inspired long ago by Erik Rosolowsky's SIGFIG:
;     http://www.cfa.harvard.edu/~erosolow/idl/lib/lib.html#SIGFIG
;
; This version written by JohnJohn Sept 29, 2005
;
; 24 Oct 2007 - If result is a single number, return scalar value
;               instead of an 1-element array. Thanks Mike Liu.
;  2 Apr 2008 - Fixed 1-element array issue, but for real this time.
;-

;;; SF_STR - The way STRING() should behave by default
function sf_str, stringin, format=format
return, strcompress(string(stringin, format=format), /rem)
end

;;; SF_TRANS_DEC - TRANSlate the DECimal point in a number of order
;;;                unity, round it, and translate back. 
function sf_trans_dec, numin, nsigin, order_inc=order_inc
nel = n_elements(numin)

;;; Double precision can't handle more than 19 sig figs
nsig = nsigin < 19

;;; Gonna have to move the decimal nsig-1 places to the right before rounding
move = nsig-1
len = max(strlen(numin))
move = move < (len-1)

;;; Create a string with just the digits, no decimal
nodec = strmid(numin, 0, 1)+strmid(numin, 2, len)

;;; Move the decimal, so nsig digits are to the left of the new
;;; decimal position
num0 = strmid(nodec,0,1+move)+'.'+strmid(nodec,1+move,len)

;;; Round the new number
num1 = sf_str(round(double(num0),/l64))
len1 = strlen(num1)

;;; If the number increases an order of magnitude after rounding, set
;;; order_inc=1 so the calling routine knows to add one to the order 
;;; of magnitude
order_inc = len1 gt nsig
;;; Move the decimal back and return to sender
num  = strmid(num1, 0, 1)+'.'+strmid(num1, 1, nsig-1)
return, num
end

function sigfig, NumIn, Nfig $
                 , string_return=string_return $
                 , scientific=scientific $
                 , numerical=numerical $
                 , plusses=plusses

Num = double(NumIn)
Nel = n_elements(Num)

;;; Convert the input number to scientific notation
TestString = sf_str(abs(double(Num)), format='(e)')
Epos = strpos(TestString[0], 'e')

;;; Test sign of the order
Osign = intarr(Nel)+1
StrOsign = strmid(TestString, Epos+1, 1)
Wneg = where(strosign eq '-', Nneg) 
if Nneg gt 0 then Osign[Wneg] = -1

;;; Test sign of numbers, form string of minus signs for negative vals
NegSign = strarr(Nel) + (keyword_set(plusses) ? '+' : '')
Negative = where(Num lt 0, Nneg)
if Nneg gt 0 then NegSign[Negative] = '-'

;;; What are the orders of magnitude of the values?
Order = fix(sf_str(strmid(TestString, Epos+2, 2)))

;;; Convert all values to order unity for rounding
NumUnit = strmid(TestString,0,epos)

;;; Use TRANS_DEC to round unit values
NumTrans = sf_trans_dec(NumUnit, Nfig, order_inc=Order_Inc)
Order = order + Osign*order_inc
Len = strlen(NumTrans[0])

;;; Exit early without looping for /NUMERICAL or /SCIENTIFIC
if keyword_set(numerical) then begin
    NumRound = NegSign+NumTrans+'e'+StrOsign+sf_str(Order)
    if n_elements(NumRound) eq 1 then return, double(NumRound[0]) else $
      return, double(NumRound)
endif
if keyword_set(scientific) then begin
    NumRound = NegSign+NumTrans+'e'+StrOsign+sf_str(Order)
    if n_elements(NumRound) eq 1 then return, NumRound[0] else $
      return, NumRound
endif

NumRound = strarr(Nel)
for i = 0, Nel-1 do begin
    if Osign[i]*Order[i]+1 gt Nfig then Format = '(I40)' else begin
        d = sf_str(fix(Nfig-(Osign[i]*Order[i])-1) > 0)
        Format = '(F40.' + d + ')'
    endelse
    New = NumTrans[i] * 10d^(Osign[i] * Order[i])
    NumRound[i] = NegSign[i]+sf_str(New, format=Format)
endfor
if n_elements(NumRound) eq 1 then return, NumRound[0]
return, NumRound
end

function strnumber, st, val, hex = hexflg
;+
; NAME:
;      STRNUMBER
; PURPOSE:
;      Function to determine if a string is a valid numeric value.
;
; CALLING SEQUENCE:
;      result = strnumber( st, [val, /HEX] )
;
; INPUTS:
;      st - any IDL scalar string
;
; OUTPUTS:
;      1 is returned as the function value if the string st has a
;      valid numeric value, otherwise, 0 is returned.
;
; OPTIONAL OUTPUT:
;      val - (optional) value of the string.  real*8
;
; OPTIONAL INPUT KEYWORD:
;       /HEX - If present and nonzero, the string is treated as a hexadecimal
;             longword integer.
;
; EXAMPLES:
;      IDL> res = strnumber(' ',val)
;           returns res=0 (not a number) and val is undefined
;
;      IDL> res = strnumber('0.2d', val)
;           returns res=1 (a valid number), and val = 0.2000d
;              
; NOTES:
;      (1) STRNUMBER was modified in February 1993 to include a special test for 
;      empty or null strings, which now returns a 0 (not a number).    Without
;      this special test, it was found that a empty string (' ') could corrupt
;      the stack.
;
;       (2) STRNUMBER will return a string such as '23.45uyrg' as a valid 
;      number (=23.45) since this is how IDL performs the type conversion.  If
;      you want a stricter definition of valid number then use the VALID_NUM
;      function.
; HISTORY:
;      version 1  By D. Lindler Aug. 1987
;      test for empty string, W. Landsman          February, 1993
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Hex keyword added.  MRG, RITSS, 15 March 2000.
;-
 if N_params() EQ 0 then begin
      print,'Syntax - result = strnumber( st, [val, /HEX] )'
      return, 0
 endif

 newstr = strtrim( st )

 if ( newstr EQ '' ) then return, 0    ;Empty string is a special case

 On_IOerror, L1                 ;Go to L1 if conversion error occurs

  If (NOT keyword_set(hexflg)) Then Begin
   val = double( newstr )
 EndIf Else Begin
   val = 0L
   reads, newstr, val, Format="(Z)"
 EndElse

 return, 1                      ;No conversion error

 L1: return, 0                  ;Conversion error occured

 end

pro STRREPLACE, Strings, Find1, Replacement1

;+
; NAME:
;        STRREPLACE
;
; PURPOSE:
;        The STRREPLACE procedure replaces the contents of one string
;        with another.  The first occurrence of the search substring, Find
;        within the source string, String is replaced by the string,
;        Replacement.
;
; CATEGORY:
;        String Processing.
;
; CALLING SEQUENCE:
;
;        STRREPLACE, String, Find, Replacement
;
; INPUTS:
;        String:   The string to have substring(s) replaced.  If String is
;                  an array, Find is replaced by Replacement in the first
;                  occurrence of Find of every element of the array.
;
;        Find:     The scalar substring to be replaced. If this argument is
;                  not a string, it is converted using IDL's default
;                  formatting rules.
;
;        Replacement:   A scalar string to replace the Find substring. If
;                  this argument is not a string, it is converted using IDL's
;                  default formattting rules.
;
; EXAMPLE:
;
;        If the variable A contains the string "IBM is fun", the
;        substring "IBM" can be replaced with the string "Microsoft"
;        by entering:
;
;        STRREPLACE, A, 'IBM', 'Microsoft'
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, June 1995.
;-

;   Check integrity of input parameter

         NP        = N_PARAMS()
         if (NP ne 3) then message,'Must be called with 3 parameters, '+$
                   'Strings, Find, Replacement'

         sz        = SIZE(Strings)
         ns        = n_elements(sz)
         if (sz(ns-2) ne 7) then message,'Parameter must be of string type.'

         Find      = STRING(Find1)
         pos       = STRPOS(Strings,Find)
         here      = WHERE(pos ne -1, nreplace)

         if (nreplace eq 0) then return

         Replacement=STRING(Replacement1)
         Flen      = strlen(Find)
         for i=0,nreplace-1 do begin

              j         = here(i)
              prefix    = STRMID(Strings(j),0,pos(j))
              suffix    = STRMID(Strings(j),pos(j)+Flen,$
                                       strlen(Strings(j))-(pos(j)+Flen))
              Strings(j) = prefix + replacement + suffix
         endfor
end

;;;****************************** readcol:
pro tap_readcol,name,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
                v16,v17,v18,v19,v20,v21,v22,v23,v24,v25, COMMENT = comment, $
                FORMAT = fmt, DEBUG=debug, SILENT=silent, SKIPLINE = skipline, $
                NUMLINE = numline, DELIMITER = delimiter
;+
; NAME:
;       READCOL
; PURPOSE:
;       Read a free-format ASCII file with columns of data into IDL vectors 
; EXPLANATION:
;       Lines of data not meeting the specified format (e.g. comments) are 
;       ignored.  Columns may be separated by commas or spaces.
;
;       Use READFMT to read a fixed-format ASCII file.   Use RDFLOAT for
;       much faster I/O (but less flexibility).    Use FORPRINT to write 
;       columns of data (inverse of READCOL).    
;
; CALLING SEQUENCE:
;       READCOL, name, v1, [ v2, v3, v4, v5, ...  v25 , COMMENT=
;           DELIMITER= ,FORMAT = , /DEBUG ,  /SILENT , SKIPLINE = , NUMLINE = ]
;
; INPUTS:
;       NAME - Name of ASCII data file, scalar string.  
;
; OPTIONAL INPUT KEYWORDS:
;       FORMAT - scalar string containing a letter specifying an IDL type
;               for each column of data to be read.  Allowed letters are 
;               A - string data, B - byte, D - double precision, F- floating 
;               point, I - integer, L - longword, Z - longword hexadecimal, 
;               and X - skip a column.
;
;               Columns without a specified format are assumed to be floating 
;               point.  Examples of valid values of FMT are
;
;       'A,B,I'        ;First column to read as a character string, then 
;                       1 column of byte data, 1 column integer data
;       'L,L,L,L'       ;Four columns will be read as longword arrays.
;       ' '             ;All columns are floating point
;
;       If a FORMAT keyword string is not supplied, then all columns are 
;       assumed to be floating point.
;
;       /SILENT - Normally, READCOL will display each line that it skips over.
;               If SILENT is set and non-zero then these messages will be 
;               suppressed.
;       /DEBUG - If this keyword is non-zero, then additional information is
;                printed as READCOL attempts to read and interpret the file.
;       COMMENT - single character specifying comment signal.   Any line 
;                beginning with this character will be skipped.   Default is
;                no comment lines.
;       DELIMITER - single character specifying delimiter used to separate 
;                columns.   Default is either a comma or a blank.
;       SKIPLINE - Scalar specifying number of lines to skip at the top of file
;               before reading.   Default is to start at the first line.
;       NUMLINE - Scalar specifying number of lines in the file to read.  
;               Default is to read the entire file
;
; OUTPUTS:
;       V1,V2,V3,...V25 - IDL vectors to contain columns of data.
;               Up to 25 columns may be read.  The type of the output vectors
;               are as specified by FORMAT.
;
; EXAMPLES:
;       Each row in a file position.dat contains a star name and 6 columns
;       of data giving an RA and Dec in sexigesimal format.   Read into IDL 
;       variables.   (NOTE: The star names must not include the delimiter 
;       as a part of the name, no spaces or commas as default.)
;
;       IDL> FMT = 'A,I,I,F,I,I,F'
;       IDL> READCOL,'position.dat',F=FMT,name,hr,min,sec,deg,dmin,dsec  
;
;       The HR,MIN,DEG, and DMIN variables will be integer vectors.
;
;       Alternatively, all except the first column could be specified as
;       floating point.
;
;       IDL> READCOL,'position.dat',F='A',name,hr,min,sec,deg,dmin,dsec 
;
;       To read just the variables HR,MIN,SEC
;       IDL> READCOL,'position.dat',F='X,I,I,F',HR,MIN,SEC
;
; RESTRICTIONS:
;       This procedure is designed for generality and not for speed.
;       If a large ASCII file is to be read repeatedly, it may be worth
;       writing a specialized reader.
;
;       Columns to be read as strings must not contain the delimiter character
;       (i.e. commas or spaces by default).   Either change the default 
;       delimiter with the DELIMITER keyword, or use READFMT to read such files.
;
;       Numeric values are converted to specified format.  For example,
;       the value 0.13 read with an 'I' format will be converted to 0.
;
; PROCEDURES CALLED
;       GETTOK(), NUMLINES(), STRNUMBER()
;
; MINIMUM IDL VERSION:
;       V5.3 (Uses STRSPLIT)
; REVISION HISTORY:
;       Written         W. Landsman                 November, 1988
;       Modified             J. Bloch                   June, 1991
;       (Fixed problem with over allocation of logical units.)
;       Added SKIPLINE and NUMLINE keywords  W. Landsman    March 92
;       Read a maximum of 25 cols.  Joan Isensee, Hughes STX Corp., 15-SEP-93.
;       Call NUMLINES() function W. Landsman          Feb. 1996
;       Added DELIMITER keyword  W. Landsman          Nov. 1999
;       Fix indexing typos (i for k) that mysteriously appeared W. L. Mar. 2000
;       Hexadecimal support added.  MRG, RITSS, 15 March 2000.
;       Default is comma or space delimiters as advertised   W.L. July 2001
;       Faster algorithm, use STRSPLIT if V5.3 or later  W.L.  May 2002
;       Accept null strings separated by delimiter ,e.g. ',,,'
;       Use SCOPE_VARFETCH instead of EXECUTE() for >V6.1  W.L. Jun 2005
;       Added compile_opt idl2   W. L.  July 2005
;-
  On_error,2                           ;Return to caller
  compile_opt idl2

  if N_params() lt 2 then begin
     print,'Syntax - READCOL, name, v1, [ v2, v3,...v25, '
     print,'        FORMAT= ,/SILENT  ,SKIPLINE =, NUMLINE = , /DEBUG]'
     return
  endif

  no_exec = !VERSION.RELEASE GE '6.1'
; Get number of lines in file

   nlines = NUMLINES( name )
   if nlines LT 0 then return

   if keyword_set(DEBUG) then $
      message,'File ' + name+' contains ' + strtrim(nlines,2) + ' lines',/INF

   if not keyword_set( SKIPLINE ) then skipline = 0
   nlines = nlines - skipline
   if keyword_set( NUMLINE) then nlines = numline < nlines

  ncol = N_params() - 1           ;Number of columns of data expected
  vv = 'v' + strtrim( indgen(ncol)+1, 2)
  nskip = 0

  if N_elements(fmt) GT 0 then begin    ;FORMAT string supplied?

    if size(fmt,/tname) NE 'STRING' then $
       message,'ERROR - Supplied FORMAT keyword must be a scalar string'
;   Remove blanks from format string
    frmt = strupcase(strcompress(fmt,/REMOVE))   
    remchar, frmt, '('                  ;Remove parenthesis from format
    remchar, frmt, ')'           

;   Determine number of columns to skip ('X' format)
    pos = strpos(frmt, 'X', 0)

    while pos NE -1 do begin
        pos = strpos( frmt, 'X', pos+1)
        nskip = nskip + 1
    endwhile

  endif else begin                     ;Read everything as floating point

    frmt = 'F'
    if ncol GT 1 then for i = 1,ncol-1 do frmt = frmt + ',F'
    if not keyword_set( SILENT ) then message, $
      'Format keyword not supplied - All columns assumed floating point',/INF

  endelse

  nfmt = ncol + nskip
  idltype = intarr(nfmt)

; Create output arrays according to specified formats

   k = 0L                                     ;Loop over output columns
   hex = bytarr(nfmt)
   for i = 0L, nfmt-1 do begin

       fmt1 = gettok( frmt, ',' )
       if fmt1 EQ '' then fmt1 = 'F'         ;Default is F format
       case strmid(fmt1,0,1) of 
          'A':  idltype[i] = 7          
          'D':  idltype[i] = 5
          'F':  idltype[i] = 4
          'I':  idltype[i] = 2
          'B':  idltype[i] = 1
          'L':  idltype[i] = 3
          'Z':  begin 
                idltype[i] = 3               ;Hexadecimal
                hex[i] = 1b
                end
          'X':  idltype[i] = 0               ;IDL type of 0 ==> to skip column
         ELSE:  message,'Illegal format ' + fmt1 + ' in field ' + strtrim(i,2)
      endcase

; Define output arrays

      if idltype[i] GT 0 then begin
          if no_exec then (SCOPE_VARFETCH(vv[k], LEVEL=0))= $
	        make_array(nlines,TYPE = idltype[i]) else $
          tst = execute(vv[k] + '= make_array(nlines,TYPE = idltype[i] )' ) 
           k = k+1
      endif

   endfor
   goodcol = where(idltype)
   idltype = idltype[goodcol]
   check_numeric = (idltype NE 7)
   openr, lun, name, /get_lun
   ngood = 0L

   temp = ' '
   if !VERSION.RELEASE GE '5.6' then skip_lun,lun,skipline, /lines else $
   if skipline GT 0 then $
       for i = 0, skipline-1 do readf, lun, temp        ;Skip any lines

   if not keyword_set(delimiter) then delimiter = ' ,'
 
   for j = 0L, nlines-1 do begin

      readf, lun, temp
      if strlen(temp) LT ncol then begin    ;Need at least 1 chr per output line
          ngood = ngood-1
          if not keyword_set(SILENT) then $
                       message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF
          goto, BADLINE 
       endif
    k = 0
    temp = strtrim(temp,1)                  ;Remove leading spaces
    if keyword_set(comment) then if strmid(temp,0,1) EQ comment then begin
          ngood = ngood-1
          if keyword_set(DEBUG) then $
                 message,'Skipping Comment Line ' + strtrim(skipline+j+1,2),/INF
          goto, BADLINE 
       endif

    var = strsplit(strcompress(temp),delimiter,/extract, /preserve_null) 
    if N_elements(var) LT nfmt then begin 
                 if not keyword_set(SILENT) then $ 
                      message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF 
                 ngood = ngood-1            
                 goto, BADLINE         ;Enough columns?
    endif
    var = var[goodcol]

    for i = 0L,ncol-1 do begin
 
           if check_numeric[i] then begin    ;Check for valid numeric data
             tst = strnumber(var[i],val,hex=hex[i])          ;Valid number?
             if tst EQ 0 then begin            ;If not, skip this line
                 if not keyword_set(SILENT) then $ 
                      message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF 
                 ngood = ngood-1
                 goto, BADLINE 
             endif
          if no_exec then $
	      (SCOPE_VARFETCH(vv[k], LEVEL=0))[ngood] = val else $
	       tst = execute(vv[k] + '[ngood] = val')

         endif else $
         if no_exec then $
	 (SCOPE_VARFETCH(vv[k], LEVEL=0))[ngood] = var[i] else $
           tst = execute(vv[k] + '[ngood] = var[i]')

      k = k+1

  endfor

BADLINE:  ngood = ngood+1

   endfor

  free_lun,lun
  if ngood EQ 0 then begin 
     message,'ERROR - No valid lines found for specified format',/INFORM
     return
  endif

  if not keyword_set(SILENT) then $
        message,strtrim(ngood,2) + ' valid lines read', /INFORM  

; Compress arrays to match actual number of valid lines
  if no_exec then begin
  for i=0,ncol-1 do $
       (SCOPE_VARFETCH(vv[i], LEVEL=0)) = $
            (SCOPE_VARFETCH(vv[i], LEVEL=0))[0:ngood-1]
 endif else begin
  for i = 0,ncol-1 do $
      tst = execute(vv[i] + '='+ vv[i]+ '[0:ngood-1]')
 endelse 

  return
  end


;;;****************************** which
pro which,proname
;Prints full filenames in IDL !path search order for a particular routine.
; proname (input string) procedure name (.pro will be appended) to find
;24-Aug-92 JAV	Create.
;10-Mar-93 JAV	Fixed bug; last directory in !path ignored; pad with ': '
on_error, 2
if n_params() lt 1 then begin
  print,'syntax: which,proname(.pro assumed)'
  retall
endif

  pathlist = '.:' + !path + ': '		;build IDL path list
  fcount = 0					;reset file counter
  il = strlen(pathlist) - 1			;length of path string
  ib = 0					;begining substring index
  ie = strpos(pathlist,':',ib)			;ending substring index
  repeat begin					;true: found path separator
    path = strmid(pathlist,ib,ie-ib)		;extract path element
    fullname = path + '/' + proname + '.pro'	;build full filename
    openr,unit,fullname,error=eno,/get_lun	;try to open file
    if eno eq 0 then begin			;true: found file
      fcount = fcount + 1			;increment file counter
      if path eq '.' then begin			;true: in current directory
	spawn,'pwd',dot				;get current working directory
	dot = dot(0)				;convert to scalar
	print,fullname + ' (. = ' + dot + ')'	;print filename + current dir
      endif else begin				;else: not in current directory
	print,fullname				;print full name
      endelse
      free_lun,unit				;close file
    endif
    ib = ie + 1					;point beyond separator
    ie = strpos(pathlist,':',ib)		;ending substring index
    if ie eq -1 then ie = il			;point at end of path string
  endrep until ie eq il				;until end of path reached
  if fcount eq 0 then begin			;true: routine not found
    print,'which: ' + proname + '.pro not found on IDL !path.'
  endif
end


;;;****************************** SHARPCORNERS
pro sharpcorners, _REF_EXTRA=_extra
;+
; NAME:
;       SHARPCORNERS
;
; PURPOSE:
;       If you look closely, especially if you make your axes thick,
;       the corners of your plots are not sharp.  (Try example below.)
;       This procedure makes sharp corners on your plot.
;
; CALLING SEQUENCE:
;       SHARPCORNERS
;
; Graphics Keywords: 
;       [, COLOR=value] [, LINESTYLE={0|1|2|3|4|5}] [, /NOCLIP] [, 
;       THICK=value]
;
; INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       Plot corners are now sharp!
;
; EXAMPLE:
;       IDL> plot, findgen(5), XTHICK=3, YTHICK=3
;       IDL> sharpcorners, THICK=3
;
; MODIFICATION HISTORY:
;   18 Oct 2002  Written by Tim Robishaw, Berkeley
;-

corners = [0,0,0,1,1,1,1,0,0,0,0,1]

for i = 0, 4 do $
  plots, [!x.window[corners[2*i  ]], !x.window[corners[2*(i+1)  ]]], $
         [!y.window[corners[2*i+1]], !y.window[corners[2*(i+1)+1]]], $
         /NORMAL, _EXTRA=_extra

end; sharpcorners

;;;****************************** HLINE
pro hline, val,_extra=extra, xlog=xlog, min=min, max=max
;+
; NAME:
;      HLINE
;     
; PURPOSE:
;      Draw a horizontal line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      HLINE, VAL
;
; INPUTS:
;
;      VAL: The y-value where the vertical line should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a horizontal line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a horizontal line at x = 4
;      IDL> plot, findgen(10)
;      IDL> hline, 4
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
;-

nv = n_elements(val)
if 1-keyword_set(min) then if keyword_set(xlog) then min = 1d-20 else $
  min = !x.crange[0]
if 1-keyword_set(max) then if keyword_set(xlog) then max = 1d20 else $
  max = !x.crange[1]
for i = 0, nv-1 do oplot,[min,max],fltarr(2)+val[i],_extra=extra 
end



;;;****************************** VLINE
pro vline, val,_extra=extra,ylog=ylog, min=min, max=max
;+
; NAME:
;      VLINE
;     
; PURPOSE:
;      Draw a vertical line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      VLINE, VAL
;
; INPUTS:
;
;      VAL: The x-value or array of x-values where the vertical
;      line(s) should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a vertical line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a vertical line at x = 0
;      IDL> plot, findgen(10)
;      IDL> vline, 5
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
;-

nv = n_elements(val)
if 1-keyword_set(min) then if keyword_set(ylog) then min = 1d-20 else $
  min = !y.crange[0]
if 1-keyword_set(max) then if keyword_set(ylog) then max = 1d20 else $
  max = !y.crange[1]

for i = 0,nv-1 do oplot,fltarr(2)+val[i],[min,max],_extra=extra 
end

;;;******************************












pro tap_loadextra
print,'Loading TAP dependencies'
end
