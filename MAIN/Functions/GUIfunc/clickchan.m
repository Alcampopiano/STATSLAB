function clickchan(gcbo,eventdata)

if  isempty(findobj(gco,'Color','g'))
    set(gco, 'Color', 'green', 'FontSize',13, 'FontWeight','bold')
    
else
    set(gco, 'Color', 'black', 'FontSize',10,'FontWeight','normal')
    
end

end




