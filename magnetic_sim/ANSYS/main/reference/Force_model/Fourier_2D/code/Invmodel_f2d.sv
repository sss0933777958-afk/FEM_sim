module invmodel_f2d(
  //=====================
  //I/O??
  //=====================
  input logic iClk100M,//????????
  input logic iRst_n,//????
  input logic iStart,//????????????????
  input logic [31:0] iFd_m [2:0],//??????fx,fy,fz
  input logic [31:0] iPosition_a [2:0],//??
  input logic [31:0] il_hat,//?????????
  input logic [31:0] iFGain_Flux,//???
  input logic [31:0] iM2A_para[8:0],//????
  input logic [31:0] iGain_ratio,//B_ghat_Phi/B_ghat_V [mV/nWb]
  output logic oFlux_avail_pls,
  output logic [31:0] oFlux [5:0],//????
  output logic [31:0] oVd [5:0]   //???? [mV]
  
);
  //=====================
  //????
  //=====================
  parameter N_STEP=10;
  logic [9:0] cnt;//???
  logic busy;//??
  logic [31:0] pbar [2:0];//????? p_bar=p/l_aht
  
  logic [31:0] sine_phi;//?sine_phi???
  logic [31:0] cosine_phi;//?cosine_phi???
  logic [31:0] cosine_theta;//?cosine_theta???
  logic [31:0] sine_theta;//?sine_theta???
  //??phi?????theta???[0 pi]
  logic [31:0] prod_M2A[2:0];//?????????
  logic [31:0] Fd_a[2:0];//acuator_frame force
  logic [31:0] Fd_za_sq;//??f_za^2
  logic [31:0] Fd_xa_ya_sq;//??sqrt(f_xa^2+f_ya^2)
  logic [31:0] F_norm ;//??sqrt(f_xa^2+f_ya^2+f_za^2)
  logic [31:0] cosine_2_time;//2cos(theta)
  logic [31:0] sine2theta;
  logic [31:0] cosine2theta;
  logic [31:0] cosine3theta;
  logic [31:0] sine3theta;
  logic [31:0] base [9:0];//??????
  logic [31:0] coef [599:0];//?????? (h*10+e)*10+k:h:??????;e:?????;k:?????
  logic [31:0] d_entry[59:0];
  logic [31:0] w [9:0];//??????x^2,y^2,z^2,2xy,2xz,2yz,2x,2y,2z,1
  logic [3:0] k;//????? 1~9
  logic [5:0] n;//?????D???? 0~59
  logic [5:0] n_a;//??????????
  logic [3:0] k_a;//????????????
  logic [5:0] n_w;
  logic [2:0] h_m, h_a;//??????????? (0~7?6?7???)
  logic [1:0] r_m, r_a;//??????????? (0~2)
  logic [31:0] lane3_sum [5:0];//???3????
  logic [31:0] Phi_hat [5:0];  //?????
  logic [31:0] scale;//???? sqrt(||f||/FGain)
  logic [31:0] Minv [35:0];//M_bar ?????? = ?*6+?
  initial $readmemh("Mbar_inv.hex", Minv, 0, 35);
  initial $readmemh("f2d_matrix.hex", coef);
  
  //====================
  //??IP??
  //====================
  logic [31:0] AddA1,AddB1,AddRst1;//???1:??1,??2,????
  logic [31:0] AddA2,AddB2,AddRst2;//???2:??1,??2,????(????)
  logic [31:0] AddA3,AddB3,AddRst3;//???3:??1,??2,????
  
  logic [31:0] SubA1,SubB1,SubRst1;//???1
  logic [31:0] SubA2,SubB2,SubRst2;//???2
  logic [31:0] SubA3,SubB3,SubRst3;//???3
  
  logic [31:0] MultA1,MultB1,MultRst1;//???1(???????)
  logic [31:0] MultA2,MultB2,MultRst2;//???2
  logic [31:0] MultA3,MultB3,MultRst3;//???3
  
  logic [31:0] DivA1,DivB1,DivRst1;//???1:???,??,?(????)
  logic [31:0] DivA2,DivB2,DivRst2;//???2:???,??,?
  
  logic [31:0] SqrtA1,SqrtRst1;//???(??17?)

  
  
  //====================
  //??IP??
  //====================
  FP_Add fp_Add1(.clk_en(1'b1),.clock(iClk100M),.dataa(AddA1),.datab(AddB1),.result(AddRst1));
  FP_Add fp_Add2(.clk_en(1'b1),.clock(iClk100M),.dataa(AddA2),.datab(AddB2),.result(AddRst2));
  FP_Add fp_Add3(.clk_en(1'b1),.clock(iClk100M),.dataa(AddA3),.datab(AddB3),.result(AddRst3));
  
  
  FP_Sub fp_Sub1(.clk_en(1'b1),.clock(iClk100M),.dataa(SubA1),.datab(SubB1),.result(SubRst1));
  FP_Sub fp_Sub2(.clk_en(1'b1),.clock(iClk100M),.dataa(SubA2),.datab(SubB2),.result(SubRst2));
  FP_Sub fp_Sub3(.clk_en(1'b1),.clock(iClk100M),.dataa(SubA3),.datab(SubB3),.result(SubRst3));
  

  FP_Mult fp_mult1 (.clk_en(1'b1),.clock(iClk100M),.dataa(MultA1),.datab(MultB1),.result(MultRst1));
  FP_Mult fp_mult2 (.clk_en(1'b1),.clock(iClk100M),.dataa(MultA2),.datab(MultB2),.result(MultRst2));
  FP_Mult fp_mult3 (.clk_en(1'b1),.clock(iClk100M),.dataa(MultA3),.datab(MultB3),.result(MultRst3));
  
  FP_Div fp_div1(.clk_en(1'b1),.clock(iClk100M),.dataa(DivA1),.datab(DivB1),.result(DivRst1));
  FP_Div fp_div2(.clk_en(1'b1),.clock(iClk100M),.dataa(DivA2),.datab(DivB2),.result(DivRst2));
  
  FP_Sqrt fp_sqrt  (.clk_en(1'b1),.clock(iClk100M),.data(SqrtA1),.result(SqrtRst1));
  
  
  //====================
  //????
  //====================
  always @(posedge iClk100M or negedge iRst_n) begin
    if (!iRst_n) begin
      oFlux_avail_pls<=0; 
      busy<=0;
      cnt<=0;
    end
    else begin
      oFlux_avail_pls<=0;  
      if (iStart) begin
        busy<=1;
        cnt<=0;
      end
    
      else if (busy) begin
        cnt<=cnt+1;
        case(cnt)
          //====================
          //???????????
          //====================
          10'd0:begin//?????????
            MultA1<=iM2A_para[0];MultB1<=iFd_m[0];
            MultA2<=iM2A_para[1];MultB2<=iFd_m[1];
            MultA3<=iM2A_para[2];MultB3<=iFd_m[2];
            DivA1<=iPosition_a[0]; DivB1<=il_hat; //x/l_hat
            DivA2<=iPosition_a[1]; DivB2<=il_hat; //y/l_hat
          end
          10'd1:begin//?????????
            MultA1<=iM2A_para[3];
            MultA2<=iM2A_para[4];
            MultA3<=iM2A_para[5];
            DivA1<=iPosition_a[2]; DivB1<=il_hat;//z/l_hat
          end
          10'd2:begin//?????????
            MultA1<=iM2A_para[6];
            MultA2<=iM2A_para[7];
            MultA3<=iM2A_para[8];
          end
          //??????????????????
          10'd6:begin
            AddA1<=MultRst1;
            AddB1<=MultRst2;
            prod_M2A[0]<=MultRst3;
          end
          10'd7:begin
            AddA2<=MultRst1;
            AddB2<=MultRst2;
            prod_M2A[1]<=MultRst3;
            pbar[0]<=DivRst1;//x_bar
            pbar[1]<=DivRst2;//y_bar
          end
          10'd8:begin
            AddA3<=MultRst1;
            AddB3<=MultRst2;
            prod_M2A[2]<=MultRst3;
            pbar[2]<=DivRst1;//z_bar
          end
          10'd9:begin//??
            MultA1<=32'h40000000; MultB1<=pbar[0];//2x
            MultA2<=32'h40000000; MultB2<=pbar[1];//2y
            MultA3<=32'h40000000; MultB3<=pbar[2];//2z
            w[9]<=32'h3F800000;                   //1
          end
          10'd10:begin//??
            MultA1<=pbar[0]; MultB1<=pbar[0];//x^2
            MultA2<=pbar[1]; MultB2<=pbar[1];//y^2
            MultA3<=pbar[2]; MultB3<=pbar[2];//z^2
          end
          
          
          //??????????
          10'd14:begin
            AddA1<=AddRst1;
            AddB1<=prod_M2A[0];
          end
          10'd15:begin
            AddA2<=AddRst2;
            AddB2<=prod_M2A[1];
            w[6]<=MultRst1; w[7]<=MultRst2; w[8]<=MultRst3;//2x,2y,2z
            MultA1<=MultRst1; MultB1<=pbar[1];//2x*y
            MultA2<=MultRst1; MultB2<=pbar[2];//2x*z
            MultA3<=MultRst2; MultB3<=pbar[2];//2y*z
          end
          10'd16:begin
            AddA3<=AddRst3;
            AddB3<=prod_M2A[2];
            w[0]<=MultRst1; w[1]<=MultRst2; w[2]<=MultRst3;//x^2,y^2,z^2
          end
          //??????
          
          10'd21:begin//??
            w[3]<=MultRst1; w[4]<=MultRst2; w[5]<=MultRst3;//2xy,2xz,2yz
          end
          
          10'd22:begin
            MultA1<=AddRst1;MultB1<=AddRst1;//???????????f_xa^2(?????6?)
            Fd_a[0]<=AddRst1;//?????????????
          end
          10'd23:begin
            MultA2<=AddRst2;MultB2<=AddRst2;//???????????f_ya^2
            Fd_a[1]<=AddRst2;//?????????????
          end
          10'd24:begin
            MultA3<=AddRst3;MultB3<=AddRst3;//???????????f_za^2
            Fd_a[2]<=AddRst3;//?????????????
          end
          10'd29:begin
            AddA1<=MultRst1;AddB1<=MultRst2;//f_xa^2+f_ya^2(???????)
          end
          
          10'd30:begin
            Fd_za_sq<=MultRst3;//f_za^2
          end
            
          10'd37:begin
            SqrtA1<=AddRst1;//sqrt(f_xa^2+f_ya^2)(?????17?)
            AddA1<=Fd_za_sq;AddB1<=AddRst1;//f_xa^2+f_ya^2+f_za^2
          end
          
          10'd45:begin
            SqrtA1<=AddRst1;//sqrt(f_xa^2+f_ya^2+f_za^2)
          end
          10'd54:begin
            DivA1<=Fd_a[0];DivB1<=SqrtRst1;//f_xa/sqrt(f_xa^2+f_ya^2)
            DivA2<=Fd_a[1];DivB2<=SqrtRst1;//f_ya/sqrt(f_xa^2+f_ya^2)
            Fd_xa_ya_sq<=SqrtRst1;
          end
          
          10'd61:begin
            cosine_phi<=DivRst1;
            sine_phi<=DivRst2;
            
          end
            
          10'd62:begin
            DivA1<=Fd_a[2];DivB1<=SqrtRst1;//f_za/sqrt(f_xa^2+f_ya^2+f_za^2)
            DivA2<=Fd_xa_ya_sq;DivB2<=SqrtRst1;//sqrt(f_xa^2+f_ya^2)/sqrt(f_xa^2+f_ya^2+f_za^2)
            F_norm<=SqrtRst1;
          end
          
          10'd63:begin//??
            DivA1<=F_norm; DivB1<=iFGain_Flux;//||f||/FGain
          end
          
            
          10'd69:begin
            if(F_norm[30:23]==8'd0)begin
              cosine_theta<=32'h3F800000;//cos(theta)=1
              sine_theta<=32'h00000000; //sin(theta)=0
              
              cosine_phi<=32'h3F800000;//cos(phi)=1
              sine_phi<=32'h00000000; //sin(phi)=0
            end
            else if (DivRst2[30:0]<31'h38D1B717)begin //sin(theta)<1e-4???+-z
              cosine_theta<=DivRst1;      sine_theta<=32'h00000000;
              cosine_phi<=32'h3F800000;   sine_phi<=32'h00000000;
              
            end
            else begin
              cosine_theta<=DivRst1;
              sine_theta<=DivRst2;
            end
          end
          
          //=========================
          //????
          //=========================
          
          10'd70:begin
            MultA1<=32'h40000000;MultB1<=cosine_theta;//2*cos(theta)
            MultA2<=cosine_theta;MultB2<=cosine_theta;//cos(theta)^2
            MultA3<=sine_theta;MultB3<=sine_theta;//sin(theta)^2
            base[0]<=32'h3F800000;//??1:1
            base[1]<=cosine_theta;//??2:cosine_theta
            SqrtA1<=DivRst1;//sqrt(||f||/FGain)
          end
          
          10'd71:begin
            MultA1<=cosine_phi;MultB1<=sine_theta;//cos_phi*sine_theta
            MultA2<=sine_phi;MultB2<=sine_theta;//sin_phi*sine_theta
            
            
          end
            
          10'd76:begin
            cosine_2_time<=MultRst1;//2*cosine_theta
            SubA1<=MultRst2;SubB1<=MultRst3;//cos2theta
            MultA1<=MultRst1;MultB1<=sine_theta;//sin2theta
          end
          
          10'd77:begin
            base[4]<=MultRst1;
            base[5]<=MultRst2;
          end
          
          10'd82:begin
            sine2theta<=MultRst1;
            
            MultA1<=MultRst1;MultB1<=cosine_2_time;//2*cos_theta*sin2theta
            MultA2<=cosine_phi;MultB2<=MultRst1;
            MultA3<=sine_phi;MultB3<=MultRst1;
          end
          
          10'd84:begin
            cosine2theta<=SubRst1;//cos2theta
            MultA2<=cosine_2_time;MultB2<=SubRst1;//2*cos_theta*cos2theta
            base[2]<=SubRst1;//???
          end
          10'd87:begin//??
            MultA3<=SqrtRst1; MultB3<=iGain_ratio;//sqrt(||f||/FGain)*???
          end
          
          10'd88:begin
            SubA1<=MultRst1;SubB1<=sine_theta;//2*cos_theta*sin2theta-sin_theta 
            base[6]<=MultRst2;
            base[7]<=MultRst3;
          end
          
          10'd93:begin//??
            scale<=MultRst3;//??????????
          end
          
          10'd90:begin
            SubA2<=MultRst2;SubB2<=cosine_theta;//2*cos_theta*cosine2theta-cosine_theta 
          end
          
          
          
          10'd96:begin
            sine3theta<=SubRst1;
            MultA1<=cosine_phi;MultB1<=SubRst1;
            MultA2<=sine_phi;MultB2<=SubRst1;
            
          end
          
          10'd98:begin
            cosine3theta<=SubRst2;
            base[3]<=SubRst2;
          end 
          //===================
          //???
          //===================
          
          10'd102:begin
            base[8]<=MultRst1;
            base[9]<=MultRst2;
            n<=6'd0;//?n?????0?????
            k<=4'd1;//?k?????1?????
            n_a<=6'd0;
            k_a<=4'd0;
            n_w<=6'd0;
            
            h_m<=3'd0; r_m<=2'd0;
            h_a<=3'd0; r_a<=2'd0;
          end
          10'd354:begin//???????????M_inv?
            h_m<=3'd0; r_m<=2'd0;
            h_a<=3'd0; r_a<=2'd0;
          end

          10'd399:begin//??? 10'd355 ?? 399
            busy<=0;
            oFlux_avail_pls<=1;
          end
         
          default: ;
        endcase
        //==============================
        //D??:???
        //==============================
      if(cnt>=10'd103 && cnt<=10'd282)begin
          MultA1<=base[k]; MultB1<=coef[n*10+k];      //D?? n
          MultA2<=base[k]; MultB2<=coef[(n+1)*10+k];  //D?? n+1
          MultA3<=base[k]; MultB3<=coef[(n+2)*10+k];  //D?? n+2
          if(n==6'd57)begin   //60?D?????
            n<=6'd0;
            k<=k+4'd1;        //??????
          end
          else begin
            n<=n+6'd3;
          end
      end
        //=====================
        //???
        //=====================
        if(cnt>=10'd109 && cnt<=10'd288)begin
          AddB1<=MultRst1;
          AddB2<=MultRst2;
          AddB3<=MultRst3;
          if(k_a==4'd0)begin
            AddA1<=coef[n_a*10];
            AddA2<=coef[(n_a+1)*10];
            AddA3<=coef[(n_a+2)*10];
         
          end 
          else begin
            AddA1<=d_entry[n_a];
            AddA2<=d_entry[n_a+1];
            AddA3<=d_entry[n_a+2];
          end
        //?????
          if(n_a==6'd57)begin
            n_a<=6'd0;
            k_a<=k_a+4'd1;
          end
          else begin
            n_a<=n_a+6'd3;
          end  
       end  
       
       if(cnt>=10'd117 && cnt<=10'd296)begin
         d_entry[n_w]<=AddRst1;
         d_entry[n_w+1]<=AddRst2;
         d_entry[n_w+2]<=AddRst3;
        if(n_w==6'd57)begin
          n_w<=6'd0;
        end
        else begin
          n_w<=n_w+6'd3;
        end
       end
       
       
       
      if(cnt>=10'd297 && cnt<=10'd318)begin
          if(h_m<3'd6)begin//6?7??????
            MultA1<=w[r_m];   MultB1<=d_entry[h_m*10+r_m];  //?1?e=0,1,2
            MultA2<=w[r_m+3]; MultB2<=d_entry[h_m*10+r_m+3];//?2?e=3,4,5
            MultA3<=w[r_m+6]; MultB3<=d_entry[h_m*10+r_m+6];//?3?e=6,7,8
          end
          if(h_m==3'd7)begin//??8???
            h_m<=3'd0;
            r_m<=r_m+2'd1;
          end
          else begin
            h_m<=h_m+3'd1;
          end
        end
        
      if(cnt>=10'd303 && cnt<=10'd324)begin
          if(h_a<3'd6)begin
            AddB1<=MultRst1; AddB2<=MultRst2; AddB3<=MultRst3;
          if(r_a==2'd0)begin//?1?????
            AddA1<=d_entry[h_a*10+9];//?1? d_e=9 ??
            AddA2<=32'h00000000;     //?2?0??
            AddA3<=32'h00000000;     //?3?0??
          end
          else begin//?2?3??8????????
            AddA1<=AddRst1; AddA2<=AddRst2; AddA3<=AddRst3;
          end
        end
        if(h_a==3'd7)begin
          h_a<=3'd0;
          r_a<=r_a+2'd1;
        end
        else begin
          h_a<=h_a+3'd1;
        end
      end
      
      
      if(cnt>=10'd327 && cnt<=10'd332)begin//?1+?2??3??
        AddA1<=AddRst1; AddB1<=AddRst2;
        lane3_sum[cnt-10'd327]<=AddRst3;
      end
      if(cnt>=10'd335 && cnt<=10'd340)begin//???3
        AddA1<=AddRst1; AddB1<=lane3_sum[cnt-10'd335];
      end
      
      
      if(cnt>=10'd343 && cnt<=10'd348)begin
        Phi_hat[cnt-10'd343]<=AddRst1;
        MultA1<=AddRst1; MultB1<=scale;//Phi_hat*????
      end
      if(cnt>=10'd349 && cnt<=10'd354)begin
        oFlux[cnt-10'd349]<=MultRst1;//????
      end
      
      //===========================
      //?M_inv??????355~368??
      //===========================
      if(cnt>=10'd355 && cnt<=10'd368)begin
        if(h_m<3'd6)begin//h_m=???(??)?6?7???
          MultA1<=Minv[h_m*6+r_m];   MultB1<=oFlux[r_m];  //?1?j=0,1
          MultA2<=Minv[h_m*6+r_m+2]; MultB2<=oFlux[r_m+2];//?2?j=2,3
          MultA3<=Minv[h_m*6+r_m+4]; MultB3<=oFlux[r_m+4];//?3?j=4,5
      end
        if(h_m==3'd7)begin
          h_m<=3'd0;
          r_m<=r_m+2'd1;
        end
        else begin
        h_m<=h_m+3'd1;
        end
      end

      //=====================
      //?M_inv??????361~374??
      //=====================
      if(cnt>=10'd361 && cnt<=10'd374)begin
        if(h_a<3'd6)begin
          AddB1<=MultRst1; AddB2<=MultRst2; AddB3<=MultRst3;
        if(r_a==2'd0)begin//?1???0??
          AddA1<=32'h00000000; AddA2<=32'h00000000; AddA3<=32'h00000000;
        end
          else begin//?2??8????????
          AddA1<=AddRst1; AddA2<=AddRst2; AddA3<=AddRst3;
        end
      end
      if(h_a==3'd7)begin
        h_a<=3'd0;
        r_a<=r_a+2'd1;
      end
      else begin
        h_a<=h_a+3'd1;
      end
      end

    //=====================
    //?????
    //=====================
    if(cnt>=10'd377 && cnt<=10'd382)begin//?1+?2??3??
      AddA1<=AddRst1; AddB1<=AddRst2;
      lane3_sum[cnt-10'd377]<=AddRst3;
    end
    if(cnt>=10'd385 && cnt<=10'd390)begin//???3
      AddA1<=AddRst1; AddB1<=lane3_sum[cnt-10'd385];
    end
    if(cnt>=10'd393 && cnt<=10'd398)begin//??Vd
      oVd[cnt-10'd393]<=AddRst1;
    end
      
      
          /*
          10'd25:begin
            
            //======================
            //??????
            //======================
            MultA1<=iFd_direction[0];MultB1<=iFd_direction[0];//r_x*r_x
            MultA2<=iFd_direction[1];MultB2<=iFd_direction[1];//r_y*r_y
            sine_phi<=iFd_direction[2];//r_z
            
            
          end
          
          10'd26:begin
            DivA1<=iPosition_a[2]; DivB1<=il_hat;//z/l_hat
          end
          
          10'd27:begin
            AddA1<=MultRst1;AddB1<=MultRst2;
          
          end
          10'd28:begin//???????
            pbar[0]<=DivRst1;//x_bar
            pbar[1]<=DivRst2;//y_bar
          end
          
          10'd29:begin
            pbar[2]<=DivRst1;//z_bar
          end
          
          10'd14:begin
            SqrtA1<=AddRst1;//sqrt(r_x^2+r_y^2)
          end
          10'd31:begin
            cosine_phi<=SqrtRst1;
            DivA1<=iFd_direction[0];DivB1<=SqrtRst1;
            DivA2<=iFd_direction[1];DivB2<=SqrtRst1;
          end
          
          10'd38:begin
            if(cosine_phi[30:23]==8'd0)begin //cos(phi)=0;???
              cosine_theta<=32'h3F800000;//cos(theta)=1???0/0??NaN
              sine_theta<=32'h00000000; //sin(theta)=0???0/0??NaN
              
            end
            else begin
              sine_phi<=DivRst2;
              cosine_phi<=DivRst1;
            end
            
          end
          default: ;
          */
        
        
      end
    end
  end
  
endmodule

