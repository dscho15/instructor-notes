%%
% FIR-filter introduktion


fa_1 = 1500
fa_2 = 2500
fs = 8000
T = 1/fs
M = 22
m = 1:M
c_m = 1./(m*pi) .* (sin(2*pi.*m*T*fa_2) - sin(2*pi.*m*T*fa_1))
c_o = 2*T*(fa_2-fa_1)
c = [flip(c_m) c_o c_m]

% koeff plot:
figure(1)
stem(abs(c))

% lav en overførselsfunktion for at plotte den med bode.
sys = tf(c, [1], T);
bode(sys)

% prøv evt. simulink

%%
% vi har Gibbs fænomener, og vi ønsker at fjerne dem / ripples i pasbåndet.
% vi kan vælge at opdatere vores designmetode med en vinduesfunktion
% c'_m = c_m * w_m, dermed a_i = c_m-i
% her har jeg valgt at bruge et hammings-vindue
% help hamming

N = M*2
n = 0:N
ham_window = 0.54 - 0.46*cos(2*pi*n/N)

% med samme ordenstal fås der:

%figure(2)
%stem(ham_window)
%tf_ham = tf(ham_window, [1], T)
%bode(tf_ham)

figure(2)
sys = tf(c.*ham_window, [1], T)
bode(sys)

%%
% obs en nærmest flad ripple
% et andet vindue er: 
% kaiser vinduet ; det kan justere side lobes og main lobes

%% 
% vi kan også approximere ordenstallet vha viden om bredden på sidelobben 
% for metoden der benyttes.
% ham har B = 2 * bredden = 2 * 2 

B = 4
df = 1000
M = B * fs / (2 * df)

% og prøve endnu en gang


