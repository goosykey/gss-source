function recipient = mailbot(recipient, message, subject, attachments)
% MATLABMAIL Send an email from a matlabbot gmail account.
%
% MATLABMAIL( recipient, message, subject, attachments )
%
%   sends the character string stored in 'message' with subjectline 'subject'
%   to the address in 'recipient' from yln.matlabbot@gmail.com.
%   'attachments' (optional) must be a string or a cell array of strings
%   containing filename(s) of the file(s) being attached.


sender = 'yln.matlabbot@gmail.com';
psswd = 'supermegapass202+BS';

if nargin < 1
    error('U MAD BRO?');
elseif nargin < 2
    message = 'Matlab bot automatic message. Have a nice day.';
elseif nargin < 3
    subject = 'MATLAB BOT >_<';    
elseif nargin < 4
    attachments = [];
end

setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(recipient, subject, message, attachments);