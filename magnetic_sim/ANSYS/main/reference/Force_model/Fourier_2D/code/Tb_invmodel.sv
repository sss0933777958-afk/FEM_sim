`timescale 1ns/1ps //xns:#????????yps:??????????#x:??x???
module tb_invmodel;
  //===========================
  //????
  //===========================
  logic iClk100M;//????????
  logic iRst_n;//????
  logic iStart;//????????????????
  logic [31:0] iFd_m [2:0] ;//??????fx,fy,fz
  logic [31:0] iPosition_a [2:0];//??
  logic [31:0] il_hat;//?????????
  logic [31:0] iFGain_Flux;//???
  logic [31:0] iM2A_para [8:0];//????
  logic [31:0] iGain_ratio;//B_ghat_Phi/B_ghat_V [mV/nWb]
  logic [31:0] oVd [5:0];//???? [mV]
  
  //===========================
  //????
  //===========================
  logic oFlux_avail_pls;
  logic [31:0] oFlux [5:0];//????
  
  //===========================
  //????
  //===========================
  initial iClk100M=0;//??????
  time t_start,t_done;//??????????
  always #5 iClk100M=~iClk100M;//?????10ns???
  parameter real    FREQ = 1000.0; //???? [Hz]
  parameter real    TS   = 10e-6;  //???? [s]???? 100 kHz
  parameter integer NS   = 300;    //??????? = 1/(FREQ*TS)

  real    t;    //???? [s]
  real    fx;   //???? [pN]
  integer k;   //??????
  integer fout; //???

  

  
  //initial ??(????????? 0 ??? begin ? end ??????????????????????????)
  initial begin
    //1.??:??0?????
    
    
    iRst_n=0;
    iStart=0;
    //????l_hat(PC???)
    iPosition_a[0]=$shortrealtobits(shortreal'(0.0));//x = 10 µm
    iPosition_a[1]=$shortrealtobits(shortreal'(0.0));//y = 20 µm
    iPosition_a[2]=$shortrealtobits(shortreal'(0.0));//z = -30 µm
    il_hat=$shortrealtobits(shortreal'(870.71870724));//l_hat µm
    iFGain_Flux=$shortrealtobits(shortreal'(2.853234e-4));//F_ghat_Phi [pN/nWb^2]    
    //?????(???)
    /*
    iFd_m[0]=$shortrealtobits(shortreal'(1.0));//Fd_xm=1pN
    iFd_m[1]=$shortrealtobits(shortreal'(0.0));//Fd_ym=1pN
    iFd_m[2]=$shortrealtobits(shortreal'(0.0));//Fd_zm=1pN
    */
    
    
    iGain_ratio=$shortrealtobits(shortreal'(20.55830));//mV/nWb
    
    
    //???????
    //xa
    iM2A_para[0]=$shortrealtobits(shortreal'(0.816497));//shortreal':?()???????shortreal??
    iM2A_para[1]=$shortrealtobits(shortreal'(0.0));
    iM2A_para[2]=$shortrealtobits(shortreal'(-0.57735));
    //ya
    iM2A_para[3]=$shortrealtobits(shortreal'(-0.408248));//shortreal':?()???????shortreal??
    iM2A_para[4]=$shortrealtobits(shortreal'(0.707107));
    iM2A_para[5]=$shortrealtobits(shortreal'(-0.577350));
    //za
    iM2A_para[6]=$shortrealtobits(shortreal'(0.408248));//shortreal':?()???????shortreal??
    iM2A_para[7]=$shortrealtobits(shortreal'(0.707107));
    iM2A_para[8]=$shortrealtobits(shortreal'(0.577350));
    
    
    repeat (3) @(posedge iClk100M);
    //2. ?? reset
    iRst_n<=1;
    repeat (2) @(posedge iClk100M);
    
    
    
    fout=$fopen("vd_sine.txt","w");//??????"w"=??



    for(k=0; k<NS; k++)begin
      t  = k*TS;                                   //?kk????
      fx = $sin(2.0*3.141592653589793*FREQ*t);      //fx = sin(2*pi*f*t)

      iFd_m[0]=$shortrealtobits(shortreal'(fx));    //??DUT
      iFd_m[1]=32'h00000000;
      iFd_m[2]=32'h00000000;

    //????? iStart??????????????
      @(posedge iClk100M); iStart<=1;   //?1?????start
      @(posedge iClk100M); iStart<=0;   //?2??????start
      repeat(998) @(posedge iClk100M);  //??998???
    
      $fwrite(fout,"%e %e %e %e %e %e %e %e\n", t, fx,
      $bitstoshortreal(oVd[0]), $bitstoshortreal(oVd[1]),
      $bitstoshortreal(oVd[2]), $bitstoshortreal(oVd[3]),
      $bitstoshortreal(oVd[4]), $bitstoshortreal(oVd[5]));
    
    end
    $fclose(fout);
    $stop;
    /*
    repeat (3) @(posedge iClk100M);//?3????
    //2. ????
    iRst_n<=1;
    repeat (2) @(posedge iClk100M);//????
    //3.??????:?1?????????0
    iStart<=1;
    @(posedge iClk100M);
    iStart <=0;
    t_start=$time;//????????
    @(posedge oFlux_avail_pls);//?oFlux_avail_pls??1
    t_done=$time;//??????
    $display("t_start = %0d ns,t_done = %0d ns",t_start,t_done);
    $display("step=%0d",(t_done-t_start)/10);
    */
    /*
    
    $display("x_bar=%f,y_bar=%f,z_bar=%f",
      $bitstoshortreal(f2d_1.pbar[0]),
      $bitstoshortreal(f2d_1.pbar[1]),
      $bitstoshortreal(f2d_1.pbar[2])
      );
    */
    //????
    $display("fd_xa=%f,fd_ya=%f,fd_za=%f",
      $bitstoshortreal(f2d_1.Fd_a[0]),
      $bitstoshortreal(f2d_1.Fd_a[1]),
      $bitstoshortreal(f2d_1.Fd_a[2])
    );
    //???????
    /*
    $display("cosine_theta=%f,sine_theta=%f,cosine_phi=%f,sine_phi=%f,
              cosine2theta=%f,cosine3theta=%f,sine2theta=%f,sine3theta=%f
              base1=%f,base2=%f,base3=%f,base4=%f,base5=%f,base6=%f,base7=%f,
              base8=%f,base9=%f,base10=%f,
              coef[0]=%f coef[599]=%e",
      $bitstoshortreal(f2d_1.cosine_theta),
      $bitstoshortreal(f2d_1.sine_theta),
      $bitstoshortreal(f2d_1.cosine_phi),
      $bitstoshortreal(f2d_1.sine_phi),
      $bitstoshortreal(f2d_1.cosine2theta),
      $bitstoshortreal(f2d_1.cosine3theta),
      $bitstoshortreal(f2d_1.sine2theta),
      $bitstoshortreal(f2d_1.sine3theta),
      
      $bitstoshortreal(f2d_1.base[0]),
      $bitstoshortreal(f2d_1.base[1]),
      $bitstoshortreal(f2d_1.base[2]),
      $bitstoshortreal(f2d_1.base[3]),
      $bitstoshortreal(f2d_1.base[4]),
      $bitstoshortreal(f2d_1.base[5]),
      $bitstoshortreal(f2d_1.base[6]),
      $bitstoshortreal(f2d_1.base[7]),
      $bitstoshortreal(f2d_1.base[8]),
      $bitstoshortreal(f2d_1.base[9]),
      
      $bitstoshortreal(f2d_1.coef[0]), 
      $bitstoshortreal(f2d_1.coef[599])
      
    );
     */
     /*
      $display("d0=%e d1=%e d2=%e d3=%e d9=%e",
        $bitstoshortreal(f2d_1.d_entry[0]), $bitstoshortreal(f2d_1.d_entry[1]),
        $bitstoshortreal(f2d_1.d_entry[2]), $bitstoshortreal(f2d_1.d_entry[3]),
        $bitstoshortreal(f2d_1.d_entry[9]));
      $display("d10=%e d30=%e d59=%e",
        $bitstoshortreal(f2d_1.d_entry[10]), $bitstoshortreal(f2d_1.d_entry[30]),
        $bitstoshortreal(f2d_1.d_entry[59]));
    */
    /*
    $display("w0=%e w1=%e w2=%e w3=%e w4=%e",
      $bitstoshortreal(f2d_1.w[0]), $bitstoshortreal(f2d_1.w[1]),
      $bitstoshortreal(f2d_1.w[2]), $bitstoshortreal(f2d_1.w[3]),
      $bitstoshortreal(f2d_1.w[4]));
    $display("w5=%e w6=%e w7=%e w8=%e w9=%e",
      $bitstoshortreal(f2d_1.w[5]), $bitstoshortreal(f2d_1.w[6]),
      $bitstoshortreal(f2d_1.w[7]), $bitstoshortreal(f2d_1.w[8]),
      $bitstoshortreal(f2d_1.w[9]));
    */
    /*
    $display("Phi1=%e Phi2=%e Phi3=%e",
      $bitstoshortreal(f2d_1.oFlux[0]), $bitstoshortreal(f2d_1.oFlux[1]),
      $bitstoshortreal(f2d_1.oFlux[2]));
    $display("Phi4=%e Phi5=%e Phi6=%e",
      $bitstoshortreal(f2d_1.oFlux[3]), $bitstoshortreal(f2d_1.oFlux[4]),
      $bitstoshortreal(f2d_1.oFlux[5]));
      */
      /*
     $display("Vd1=%e Vd2=%e Vd3=%e",
  $bitstoshortreal(f2d_1.oVd[0]), $bitstoshortreal(f2d_1.oVd[1]),
  $bitstoshortreal(f2d_1.oVd[2]));
  $display("Vd4=%e Vd5=%e Vd6=%e",
  $bitstoshortreal(f2d_1.oVd[3]), $bitstoshortreal(f2d_1.oVd[4]),
  $bitstoshortreal(f2d_1.oVd[5])); 
      
    $stop;
    */
  end
  //=======================
  //????
  //=======================
  invmodel_f2d f2d_1(

    .iStart (iStart),
    .iClk100M (iClk100M),
    .iRst_n (iRst_n),
    .iFd_m (iFd_m),
    .iM2A_para (iM2A_para),
    .iPosition_a (iPosition_a),
    .il_hat (il_hat),
    .iFGain_Flux (iFGain_Flux),
    .oFlux_avail_pls (oFlux_avail_pls),
    .oFlux (oFlux),
    .iGain_ratio(iGain_ratio),
    .oVd (oVd)
  );
endmodule