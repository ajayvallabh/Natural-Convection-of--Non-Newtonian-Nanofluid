function MakeVideo()

 %       current_dir=pwd;
 %       addpath( strcat(current_dir,'/matlab/basics'), strcat(current_dir,'/matlab/tikz'));
        close all ;

        Ra = 170822.0; 
        Pr = 696555.5 ; 
        values = strcat('Ra',num2str(Ra,'%.1f'),'_Pr',num2str(Pr,'%.1f')) ; 
        P = load( strcat('F:\Ajay\maincode/output_u_theta_psi_',values,'.txt' ) ) ;
		
        RenderVideo = 1 ; 
        
        if RenderVideo
            Video = VideoWriter(strcat('F:/video_psi_theta_',values,'.avi'),'Uncompressed AVI' );   
            Video.FrameRate = 10;
            open(Video);
        end
				
        delta_x = 0.01 ;
        n = ceil( 1/delta_x ) + 1   ; % number of elements + 1 = number of nodes. 
        m = 1 * (n-1) + 1 ; 
        X_oneCol = P(1,1:m ) ;
        Y_oneCol = P(2,1:m:end ) ;
        [X,Y] = meshgrid( X_oneCol, Y_oneCol ); 
%-----------------------------------------------------------------------        

%       Plot streamfunction contour (streamlines)         

        if RenderVideo
            FIG = figure('visible','off');
        else
            FIG = figure() ; 
        end        
        
        subplot(3,1,1); hold on; grid off; box on; %axis equal; 
        xlabel('\textbf{X}','interpreter','latex');
        ylabel('\textbf{Y}','interpreter','latex');
        xlim([0,1]); 
        ylim([0,1]);   
        daspect auto; 
        
        values = strcat('Ra=',num2str(Ra,'%.1f'),', Pr=',num2str(Pr,'%.1f')) ; 

        title( strcat( 'Stream function for ', {' '} , values ) );
        
        Z = griddata( P(1,:), P(2,:), P(5,:), X, Y ) ;  
%         caxis([0,0.003]) ;        
%         caxis([-3,3]) ;        
        [C,h_psi]= contour( X, Y, Z , 'LineColor', 'none', 'showtext', 'off', 'fill', 'on');
        colorbar;
%         v = 0:0.001:0.005 ; 
%         clabel(C,h_psi,v,'FontSize',8);
%         clabel(C,h_psi,'FontSize',8);
        
% %---------------------------------------------------------------------
%       Plot velocity profiles at x = 1, 7, 15 
        
%         x = [10] ;
%         num_x = numel(x) ; 
%         h1 = cell(1, num_x) ; 
%         h2 = cell(1, num_x) ;
%         
%         for i = 1:num_x 
%             U = P(3, (x(i)*n):m:end ) ; 
%             [ h1{i}, h2{i} ] = plot_velocity_profile(U, x(i), Y_oneCol);
%         end

% %-----------------------------------------------------------------------        

%       Plot Temperature profile 

        subplot(3,1,2) ; hold on; grid off; box on; %axis equal; 
        xlabel('\textbf{X}','interpreter','latex');
        ylabel('\textbf{Y}','interpreter','latex');
        xlim([0,1]); 
        ylim([0,1]);  
        daspect auto;
        title( strcat( 'Temperature for ', {' '} , values) );
        
        Z = griddata( P(1,:), P(2,:), P(4,:), X, Y ) ;  
        cmap = colormap(gray) ;            
%         colormap( cmap(end:-1:1,:) ) ;
%         caxis([0,1]) ;
        [~,h_theta]=contour( X, Y, Z, 'LineColor', 'none', 'showtext', 'off', 'fill', 'on');    
        colorbar;  

%-----------------------------------------------------------------------        

%       Plot Vorticity 

        subplot(3,1,3); hold on; grid off; box on; %axis equal; 
        xlabel('\textbf{X}','interpreter','latex');
        ylabel('\textbf{Y}','interpreter','latex');
        xlim([0,1]); 
        ylim([0,1]);  
        daspect auto;
        title( strcat( 'Vorticity for ', {' '} , values ) );
        daspect auto;
        Z = griddata( P(1,:), P(2,:), P(7,:), X, Y ) ;  
        cmap = colormap(gray) ;            
%         colormap( cmap(end:-1:1,:) ) ;
        
%          caxis([-0.12,0.07]) ;
%          caxis([-100,100]) ;
        [~,h_vorticity]=contour( X, Y, Z, 'LineColor', 'none', 'showtext', 'off', 'fill', 'on');    
        
        colorbar; 
        
% --------------------------------------------        
% Animation        
        
        if RenderVideo
            for j = 4:5:size(P,1)             
                Z = griddata( P(1,:), P(2,:), P(j,:), X, Y ) ;  
                set(h_theta, 'zdata', Z ) ;
                
                Z = griddata( P(1,:), P(2,:), P(j+1,:), X, Y ) ;  
                set(h_psi, 'zdata', Z ) ; 
                
                Z = griddata( P(1,:), P(2,:), P(j+3,:), X, Y ) ;  
                set(h_vorticity, 'zdata', abs(Z) ) ; 
                
%                 subplot(3,1,1) ;
%                 for i = 1: num_x 
%                     U = P(j-1, (x(i)*n):m:end ) ;                 
%                     delete( h1{i} ) ;
%                     delete( h2{i} ) ; 
%                     [ h1{i}, h2{i} ] = plot_velocity_profile(U, x(i), Y_oneCol);                
%                 end    

                cdata = print('-r300','-RGBImage') ;
                writeVideo(Video,cdata);  
                if ~rem(j, 100)
                    fprintf('%d \n', j) ;
                end
    % 			pause(0.01) ;
            end
            close(Video); 
        else
            for j = 4:5:size(P,1)             
                Z = griddata( P(1,:), P(2,:), P(j,:), X, Y ) ;  
                set(h_theta, 'zdata', Z ) ;
                
                Z = griddata( P(1,:), P(2,:), P(j+1,:), X, Y ) ;  
                set(h_psi, 'zdata', Z ) ; 
                
                Z = griddata( P(1,:), P(2,:), P(j+3,:), X, Y ) ;  
                set(h_vorticity, 'zdata', Z ) ; 

%                 subplot(3,1,1) ;
%                 for i = 1: num_x 
%                     U = P(j-1, (x(i)*n):m:end ) ;                 
%                     delete( h1{i} ) ;
%                     delete( h2{i} ) ; 
%                     [ h1{i}, h2{i} ] = plot_velocity_profile(U, x(i), Y_oneCol);                
%                 end    
             
    			pause(0.01) ;
            end            
        end
         
end

function [h1, h2] = plot_velocity_profile( U, x, Y )
    COLOR = 0.5*[1,1,1] ; 
    h1 = plot( [ x*ones(size(U)), U + x ], [ Y, Y ], 'color', COLOR ) ;
    U = U(1:4:end);
    Y = Y(1:4:end);
    h2 = plot( [ x*ones(size(U)); U + x ], [ Y; Y ], 'color', COLOR ) ;
end

