function J = TAC(x)
global br_p bz_fr deb_p

br_p = x(1);
deb_p = x(2);
bz_fr = x(3);

%linking aspen to matlab
Aspen=actxserver('Apwn.Document.40.0');
[stat,mess]=fileattrib;
Aspen.invoke('InitFromArchive2', [mess.Name '\L10_Latest.bkp']);
Aspen.Visible = 0;
Aspen.SuppressDialogs = 1;

%taking value from aspen
Aspen.Application.Tree.FindNode('\Data\Blocks\C2\Subobjects\Design Specs\1\Input\VALUE\1').value = br_p;
Aspen.Application.Tree.FindNode('\Data\Blocks\C3\Subobjects\Design Specs\1\Input\VALUE\1').value = deb_p;
Aspen.Application.Tree.FindNode("\Data\Blocks\MKP\Input\OUTMLFL").value = bz_fr;

%running simulation
Run2(Aspen.Engine);

%catch error
status = Aspen.Application.Tree.FindNode('\Data\Results Summary\Run-Status\Output\PER_ERROR').value;

if status==0
    %Reb duty of 3 columns
    C1Q = Aspen.Application.Tree.FindNode("\Data\Blocks\C1\Output\REB_DUTY").value;
    C2Q = Aspen.Application.Tree.FindNode("\Data\Blocks\C2\Output\REB_DUTY").value;
    C3Q = Aspen.Application.Tree.FindNode("\Data\Blocks\C3\Output\REB_DUTY").value;
    hot = (C2Q)*0.0036*24*300*8.22 + (C1Q)*0.0036*24*300*9.8 + (C3Q)*0.0036*24*300*9.8;

    %Condensor duty of 3 columns
    C1_CQ = Aspen.Application.Tree.FindNode("\Data\Blocks\C1\Output\COND_DUTY").value;
    C2_CQ = Aspen.Application.Tree.FindNode("\Data\Blocks\C2\Output\COND_DUTY").value;
    C3_CQ = Aspen.Application.Tree.FindNode("\Data\Blocks\C3\Output\COND_DUTY").value;
    cold = (C1_CQ+C2_CQ+C3_CQ)*0.0036*24*300*0.354;
    
    %Operating cost
    opx = hot - cold;

    %no of stages, dia, temp of top & bottom of 3 columns
    C1N = Aspen.Application.Tree.FindNode("\Data\Blocks\C1\Input\NSTAGE").value;
    C1N = double(C1N);
    C1_dia = Aspen.Application.Tree.FindNode("\Data\Blocks\C1\Subobjects\Column Internals\INT-1\Output\CA_DIAM6\INT-1\CS-2").value;
    
    C1_TOPT = Aspen.Application.Tree.FindNode("\Data\Streams\C1D\Output\TEMP_OUT\MIXED").value;
    C1_BOTT = Aspen.Application.Tree.FindNode("\Data\Streams\C1B\Output\TEMP_OUT\MIXED").value;
    
    C2N = Aspen.Application.Tree.FindNode("\Data\Blocks\C2\Input\NSTAGE").value;
    C2N = double(C2N);
    C2_dia = Aspen.Application.Tree.FindNode("\Data\Blocks\C2\Subobjects\Column Internals\INT-1\Output\CA_DIAM6\INT-1\CS-2").value;
    
    C2_TOPT = Aspen.Application.Tree.FindNode("\Data\Streams\C2D\Output\TEMP_OUT\MIXED").value;
    C2_BOTT = Aspen.Application.Tree.FindNode("\Data\Streams\C2B\Output\TEMP_OUT\MIXED").value;
    
    C3N = Aspen.Application.Tree.FindNode("\Data\Blocks\C3\Input\NSTAGE").value;
    C3N = double(C3N);
    C3_dia = Aspen.Application.Tree.FindNode("\Data\Blocks\C3\Subobjects\Column Internals\INT-1\Output\CA_DIAM6\INT-1\CS-2").value;
    
    C3_TOPT = Aspen.Application.Tree.FindNode("\Data\Streams\C3D\Output\TEMP_OUT\MIXED").value;
    C3_BOTT = Aspen.Application.Tree.FindNode("\Data\Streams\C3B\Output\TEMP_OUT\MIXED").value;
    
    %shape factor (column 1)
    Fq1 = 10^(0.477+0.085*log(C1N)-0.347*log(C1N)^2);
    %estimating tray cost (column 1)
    A1 = pi*(C1_dia^2)/4;
    CP1 = 10^(2.994+0.446*log(A1)+0.396*log(A1)^2);
    TC1 = CP1*C1N*1.8*Fq1;
    
    %shape factor (column 2)
    Fq2 = 1;
    %estimating tray cost (column 2)
    A2 = pi*(C2_dia^2)/4;
    CP2 = 10^(2.994+0.446*log(A2)+0.396*log(A2)^2);
    TC2 = CP2*C2N*1.8*Fq2;
    
    %shape factor (column 3)
    Fq3 = 1;
    %estimating tray cost (column 3)
    A3 = pi*(C3_dia^3)/4;
    CP3 = 10^(2.994+0.446*log(A3)+0.396*log(A3)^2);
    TC3 = CP3*C3N*1.8*Fq3;
    
    %Total Capital Cost
    %Column cost + Trays cost + Reboiler cost + Condensor cost
    C1_CC = TC1 + 17640*C1_dia^1.066*((C1N-2)*2*0.3048*1.2)^0.802+7296*(-C1_CQ/(0.568*(C1_TOPT-298)))^0.65+7296*(C1Q/(0.852*(527-C1_BOTT)))^0.65;
    C2_CC = TC2 + 17640*C2_dia^1.066*((C2N-2)*2*0.3048*1.2)^0.802+7296*(-C2_CQ/(0.568*(C2_TOPT-298)))^0.65+7296*(C2Q/(0.852*(457-C2_BOTT)))^0.65;
    C3_CC = TC3 + 17640*C3_dia^1.066*((C3N-2)*2*0.3048*1.2)^0.802+7296*(-C3_CQ/(0.568*(C3_TOPT-298)))^0.65+7296*(C3Q/(0.852*(527-C3_BOTT)))^0.65;
    
    %estimating cost of reactor
    vol_CSTR = 200;
    D_CSTR = (2*vol_CSTR/pi)^0.333;
    L_CSTR = 2*D_CSTR;
    CSTR_C = 17460*(D_CSTR^1.066)*(L_CSTR^0.802);
    cpx = C1_CC + C2_CC+ C3_CC + 2*CSTR_C;
    
    J = cpx/3 + opx;

    else

        J = 10e10;
    end
end