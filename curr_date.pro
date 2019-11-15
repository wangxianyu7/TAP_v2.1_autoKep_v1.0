function curr_date, format=format
  
  now = systime(0)
  now_struct = {year: 0 $
                , month_num: 0 $
                , day: 0 $
                , hh: 0, mm: 0, ss: 0 $
                , month_txt: 'a', weekday: 'a' }
  
  now_arr = strsplit(now,' ',/extract)
  now_time = strsplit(now_arr[3],':',/extract)
  
  now_struct.year = now_arr[4]
  now_struct.month_txt = now_arr[1]
  now_struct.day = now_arr[2]
  now_struct.mm = now_time[1]
  now_struct.ss = now_time[2]
  now_struct.hh = now_time[0]

  case now_struct.month_txt of
     'Jan' : now_struct.month_num = '1'
     'Feb' : now_struct.month_num = '2'
     'Mar' : now_struct.month_num = '3'
     'Apr' : now_struct.month_num = '4'
     'May' : now_struct.month_num = '5'
     'Jun' : now_struct.month_num = '6'
     'Jul' : now_struct.month_num = '7'
     'Aug' : now_struct.month_num = '8'
     'Sep' : now_struct.month_num = '9'
     'Oct' : now_struct.month_num = '10'
     'Nov' : now_struct.month_num = '11'
     'Dec' : now_struct.month_num = '12'
  endcase
  
  ;; print_struct,now_struct
  
  if keyword_set(format) then begin
     case format of
        'yyyy.mmdd'  :  time = string(now_struct.year,format='(i4.4)') + '.' $
                               + string(now_struct.month_num,format='(i2.2)') $
                               + string(now_struct.day,format='(i2.2)')
        'yyyymmdd'   :  time = string(now_struct.year,format='(i4.4)')  $
                               + string(now_struct.month_num,format='(i2.2)') $
                               + string(now_struct.day,format='(i2.2)')   
        'hh:mm:ss yyyymmdd'   :  time = string(now_struct.hh,format='(i2.2)') $
                                        + ":" + string(now_struct.mm,format='(i2.2)') $
                                        + ":" + string(now_struct.ss,format='(i2.2)') $
                                        + " " + string(now_struct.year,format='(i4.4)')  $
                                        + string(now_struct.month_num,format='(i2.2)') $
                                        + string(now_struct.day,format='(i2.2)')   
        'hh:mm:ss'   :  time = string(now_struct.hh,format='(i2.2)') $
                                        + ":" + string(now_struct.mm,format='(i2.2)') $
                                        + ":" + string(now_struct.ss,format='(i2.2)') 
        'yyyymmdd_hhmm'   :  time = string(now_struct.year,format='(i4.4)')  $
                                    + string(now_struct.month_num,format='(i2.2)') $
                                    + string(now_struct.day,format='(i2.2)') $
                                    + "_" + string(now_struct.hh,format='(i2.2)') $
                                    + string(now_struct.mm,format='(i2.2)') 
        'ELSE'       : print, 'curr_date '+format+ ' format not supported'
     ENDCASE
     return,time
  endif
  
  return, now_struct 
end






