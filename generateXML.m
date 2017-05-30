function generateXML(p,xtraj,utraj,K)

stateName = p.getStateFrame.getCoordinateNames;
inputName = p.getInputFrame.getCoordinateNames;
ts = xtraj.getBreaks();

docNode = com.mathworks.xml.XMLUtils.createDocument('drake');

xtraj_node = docNode.createElement('xtraj');
docNode.getDocumentElement.appendChild(xtraj_node);

entry_node = docNode.createElement(['t','0']);
docNode.getDocumentElement.appendChild(entry_node);
num_node = docNode.createElement('num');
num_text = docNode.createTextNode(num2str(length(ts)));
num_node.appendChild(num_text);
entry_node.appendChild(num_node);
xtraj_node.appendChild(entry_node);

for z = 1:length(ts)-1
    entry_node = docNode.createElement(['t',num2str(z)]);
    docNode.getDocumentElement.appendChild(entry_node);
    time_node = docNode.createElement('t');
    time_text = docNode.createTextNode(num2str(ts(z)));
    time_node.appendChild(time_text);
    entry_node.appendChild(time_node);
    
    for i=1:9
        sn = stateName{i};
        for j=1:6
            trajName = [sn,num2str(j)];
            thisElement = docNode.createElement(trajName);
            thisElement.appendChild(docNode.createTextNode(sprintf('%d',xtraj.pp.coefs(9*(z-1)+i,j))));
            entry_node.appendChild(thisElement);            
        end
    end
    xtraj_node.appendChild(entry_node);
end

utraj_node = docNode.createElement('utraj');
docNode.getDocumentElement.appendChild(utraj_node);

entry_node = docNode.createElement(['t','0']);
docNode.getDocumentElement.appendChild(entry_node);
num_node = docNode.createElement('num');
num_text = docNode.createTextNode(num2str(length(ts)));
num_node.appendChild(num_text);
entry_node.appendChild(num_node);
utraj_node.appendChild(entry_node);

for z = 1:length(ts)-1
    entry_node = docNode.createElement(['t',num2str(z)]);
    docNode.getDocumentElement.appendChild(entry_node);
    time_node = docNode.createElement('t');
    time_text = docNode.createTextNode(num2str(ts(z)));
    time_node.appendChild(time_text);
    entry_node.appendChild(time_node);
    
    for i=1:4
        sn = inputName{i};
        for j=1:2
            trajName = [sn,num2str(j)];
            thisElement = docNode.createElement(trajName);
            thisElement.appendChild(docNode.createTextNode(sprintf('%d',utraj.pp.coefs(4*(z-1)+i,j))));
            entry_node.appendChild(thisElement);            
        end
    end
    utraj_node.appendChild(entry_node);
end

ts = K.D.getBreaks;
K_node = docNode.createElement('K');
docNode.getDocumentElement.appendChild(K_node);

entry_node = docNode.createElement(['t','0']);
docNode.getDocumentElement.appendChild(entry_node);
num_node = docNode.createElement('num');
num_text = docNode.createTextNode(num2str(length(ts)));
num_node.appendChild(num_text);
entry_node.appendChild(num_node);
K_node.appendChild(entry_node);

for z = 1:length(ts)-1
    entry_node = docNode.createElement(['t',num2str(z)]);
    docNode.getDocumentElement.appendChild(entry_node);
    time_node = docNode.createElement('t');
    time_text = docNode.createTextNode(num2str(ts(z)));
    time_node.appendChild(time_text);
    entry_node.appendChild(time_node);

    for i=1:9
        for j=1:4
            KName = ['K',num2str(j),num2str(i),];
            KijElement = docNode.createElement(KName);
            for k=1:10
                K_coefs_name = ['K',num2str(k)];
                thisElement = docNode.createElement(K_coefs_name);
                thisElement.appendChild(docNode.createTextNode(sprintf('%d',K.D.pp.coefs(36*(z-1)+4*(i-1)+j,k))));
                KijElement.appendChild(thisElement);            
            end
            entry_node.appendChild(KijElement);
        end
    end
    K_node.appendChild(entry_node);
end

xmlFileName = ['aaa','.xml'];
xmlwrite(xmlFileName,docNode);
% type(xmlFileName);

end