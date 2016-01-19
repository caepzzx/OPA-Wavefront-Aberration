tic
 for j=1:num
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);  
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
       
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
        
             v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);  
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
       
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
 end
 

toc

tic

 parfor j=1:num
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);  
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
       
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
             v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);  
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
        
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
       
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,-h/2,P_w,S_w(j),I_w(j),K_con_wav(K_con,j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:); 
 end

toc