
while (1)
dmsimport('PK26_C7_1-PerpHys.VHD');
dmsimport('PK26_C7_1-PerpHys_2.VHD');
%PK26_C7_1_PerpHys_2.Minor = 1;

%PK26_C7_1_PerpHys.Minor = 1;
dms;
plot(PK26_C7_1_PerpHys.Field,PK26_C7_1_PerpHys.RawCorrX,PK26_C7_1_PerpHys_2.Field,PK26_C7_1_PerpHys_2.RawCorrX,'.-','DisplayName','PK26_C7_1_PerpHys.RawCorrX vs. PK26_C7_1_PerpHys.Field','XDataSource','PK26_C7_1_PerpHys.Field','YDataSource','PK26_C7_1_PerpHys.RawCorrX');figure(gcf)
hold all;
plot(PK26_C7_1_PerpHys_2.Analysis.Field,PK26_C7_1_PerpHys_2.Analysis.RawHcLine2,'DisplayName','PK26_C7_1_PerpHys_2.Analysis.RawHcLine2','XDataSource','PK26_C7_1_PerpHys_2.Analysis.Field','YDataSource','PK26_C7_1_PerpHys_2.Analysis.RawHcLine2');hold all;plot(PK26_C7_1_PerpHys_2.Analysis.Field,PK26_C7_1_PerpHys_2.Analysis.RawHcLine1,'DisplayName','PK26_C7_1_PerpHys_2.Analysis.RawHcLine1','XDataSource','PK26_C7_1_PerpHys_2.Analysis.Field','YDataSource','PK26_C7_1_PerpHys_2.Analysis.RawHcLine1');plot(PK26_C7_1_PerpHys_2.Analysis.Field,PK26_C7_1_PerpHys_2.Analysis.RawMsLine1,'DisplayName','PK26_C7_1_PerpHys_2.Analysis.RawMsLine1','XDataSource','PK26_C7_1_PerpHys_2.Analysis.Field','YDataSource','PK26_C7_1_PerpHys_2.Analysis.RawMsLine1');plot(PK26_C7_1_PerpHys_2.Analysis.Field,PK26_C7_1_PerpHys_2.Analysis.RawMsLine2,'DisplayName','PK26_C7_1_PerpHys_2.Analysis.RawMsLine2','XDataSource','PK26_C7_1_PerpHys_2.Analysis.Field','YDataSource','PK26_C7_1_PerpHys_2.Analysis.RawMsLine2');hold off;figure(gcf);
axis([-11000 11000 -100 100]);
hold off;
pause(5);
end