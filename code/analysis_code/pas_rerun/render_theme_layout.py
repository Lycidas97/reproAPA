"""Panel-boundary layout and artist-level receipts for the locked supplement grid."""
import json,copy
import numpy as np
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import PathCollection
from render_theme_studies import convert_text,DATASETS,swatch_string

GAP=4.0
INTERNAL_GAP=2.8
LETTER=4.4
PAD=.45
CONTENT_PAD=1.15  # Includes the invisible ascender/descender margin of fixed-size text.

def extent(r,artist,renderer):
    b=artist.get_window_extent(renderer).transformed(r.fig.dpi_scale_trans.inverted())
    return np.array([b.x0*25.4,r.height-b.y1*25.4,b.x1*25.4,r.height-b.y0*25.4])

def axisbox(r,a):
    p=a.get_position();return np.array([p.x0*160,(1-p.y1)*r.height,p.x1*160,(1-p.y0)*r.height])

def setbox(r,a,b):
    a.set_position([b[0]/160,1-b[3]/r.height,(b[2]-b[0])/160,(b[3]-b[1])/r.height])

def axis_contentbox(r,a,renderer):
    """Visible content, including annotations; an axis-off canvas is not content."""
    if a.axison:
        b=a.get_tightbbox(renderer).transformed(r.fig.dpi_scale_trans.inverted())
        return np.array([b.x0*25.4,r.height-b.y1*25.4,b.x1*25.4,r.height-b.y0*25.4])
    boxes=[]
    for artist in list(a.texts)+list(a.lines)+list(a.patches)+list(a.collections)+list(a.images):
        if not artist.get_visible():continue
        if isinstance(artist,Text) and not artist.get_text():continue
        b=artist.get_tightbbox(renderer)
        if b is None or not np.isfinite(b.extents).all():
            if not hasattr(artist,'get_datalim'):continue
            b=artist.get_datalim(a.transData).transformed(a.transData)
        if not np.isfinite(b.extents).all():continue
        b=b.transformed(r.fig.dpi_scale_trans.inverted())
        boxes.append([b.x0*25.4,r.height-b.y1*25.4,b.x1*25.4,r.height-b.y0*25.4])
    if not boxes:return axisbox(r,a)
    return np.array([min(b[0] for b in boxes),min(b[1] for b in boxes),max(b[2] for b in boxes),max(b[3] for b in boxes)])

def unionbox(boxes):
    return np.array([min(b[0] for b in boxes),min(b[1] for b in boxes),max(b[2] for b in boxes),max(b[3] for b in boxes)])

def pack(rects,target,gap=GAP):
    """Keep the slicing topology and relative allocations; use one physical gutter."""
    if len(rects)==1:return {next(iter(rects)):np.array(target,float)}
    for d in [1,0]:
        starts=sorted(set(v[d] for v in rects.values()))
        for start in starts[1:]:
            before={k:v for k,v in rects.items() if v[d+2]<=start+.01}
            after={k:v for k,v in rects.items() if v[d]>=start-.01}
            if len(before)+len(after)!=len(rects) or not before or not after:continue
            low=min(v[d] for v in before.values());high=max(v[d+2] for v in before.values());lo2=min(v[d] for v in after.values());hi2=max(v[d+2] for v in after.values())
            frac=(high-low)/(high-low+hi2-lo2);cut=target[d]+(target[d+2]-target[d]-gap)*frac
            a=np.array(target,float);c=a.copy();a[d+2]=cut;c[d]=cut+gap
            return dict(pack(before,a,gap),**pack(after,c,gap))
    raise ValueError(('Grid is not a slicing layout',rects))

def signature(h,label):
    parts=[label]
    for name in ['get_color','get_linestyle','get_marker','get_facecolor','get_edgecolor']:
        if hasattr(h,name):parts.append(str(getattr(h,name)()))
    return tuple(parts)

def add_key(r,axes,rect,top_offset=0):
    """Keys belong to one final panel; shared-y components share only that key."""
    entries=[];seen=set();title=''
    for ax in axes:
        lg=ax.get_legend()
        if lg is None:continue
        if lg.get_title().get_text():title=lg.get_title().get_text()
        for h,t in zip(lg.legend_handles,lg.get_texts()):
            label=getattr(t,'dataset_original',t.get_text());sig=signature(h,label)
            if sig not in seen:entries.append((h,label));seen.add(sig)
        lg.remove()
    names=[k for k in DATASETS if any(k in r.dataset_axis_names.get(ax,[]) for ax in axes)]
    legends=[];yy=rect[1]+PAD+top_offset;available=rect[2]-rect[0]-LETTER-PAD*2
    for key_entries,fontsize,dataset_key in [(entries,7,False),([(Patch(facecolor=DATASETS[k][0],edgecolor='none'),k) for k in names],8,True)]:
        if not key_entries:continue
        rows=[[]];width=0
        for h,label in key_entries:
            text=label if dataset_key else swatch_string(label)[0]
            if dataset_key and r.mmwidth(text,fontsize)+7>available:
                text=text.replace('–','–\n',1)
            ww=r.mmwidth(text,fontsize)+7.0
            if width+ww>available and rows[-1]:rows.append([]);width=0
            rows[-1].append((h,text if dataset_key else label));width+=ww
        for row in rows:
            handles,labels=zip(*row)
            lg=r.fig.legend(handles,labels,loc='upper left',bbox_to_anchor=((rect[0]+LETTER+PAD)/160,1-yy/r.height),bbox_transform=r.fig.transFigure,ncol=len(row),frameon=False,fontsize=fontsize,handlelength=1.1 if dataset_key else 1.65,handletextpad=.4,columnspacing=.7,labelspacing=.2,borderpad=0,borderaxespad=0,title=title if not legends and not dataset_key else None,title_fontsize=7)
            if not dataset_key:
                for t in lg.get_texts():convert_text(t)
            lg._dataset_key=dataset_key
            r.fig.canvas.draw();hh=extent(r,lg,r.fig.canvas.get_renderer())[3]-yy
            legends.append(lg);yy+=hh+.65
    return legends,yy-rect[1]+.45 if legends or top_offset else 0,names

def text_items(ax):
    ts=list(ax.texts)+[ax.title,ax.xaxis.label,ax.yaxis.label]
    if ax.axison:
        for values,labels,lim in [(ax.get_xticks(),ax.get_xticklabels(),ax.get_xlim()),(ax.get_yticks(),ax.get_yticklabels(),ax.get_ylim())]:
            ts.extend(t for v,t in zip(values,labels) if min(lim)<=v<=max(lim))
        ts += [ax.xaxis.get_offset_text(),ax.yaxis.get_offset_text()]
    return [t for t in ts if t.get_visible() and t.get_text()]

def layout(r,groups):
    old={p['panel']:np.array([p['rectangle_mm'][0],p['rectangle_mm'][1],p['rectangle_mm'][0]+p['rectangle_mm'][2],p['rectangle_mm'][1]+p['rectangle_mm'][3]]) for p,aa in groups}
    root=np.array([1.5,1.5,158.5,max(v[3] for v in old.values())])
    if r.compact_plan:
        requested={r.fid+k:np.array([x,y,x+w,y+h]) for k,(x,y,w,h) in r.compact_plan.items()}
        assert set(requested)==set(old), (r.fid,'Compact plan must place every existing final panel exactly once')
        ordered=sorted(groups,key=lambda pair:(requested[pair[0]['panel']][1],requested[pair[0]['panel']][0]))
        new={};renamed={}
        for i,(p,axes) in enumerate(ordered):
            previous=p['panel'];p['previous_panel']=previous;p['panel']=r.fid+chr(65+i)
            if 'redraw-v13-' in r.out.name or ('redraw-v14-' in r.out.name and r.fid!='D03'):p['previous_panel']=p['panel']
            new[p['panel']]=requested[previous];renamed[p['panel']]=old[previous]
        groups=ordered;old=renamed
    else:new=pack(old,root)
    r.strict_groups=[];r.strict_text_artists=[];r.panel_keys=[];r.content_headers={};r.legend_ownership=[]
    for p,axes in groups:
        axes=sorted(axes,key=lambda a:r.fig.axes.index(a))
        rect=new[p['panel']];oldrect=old[p['panel']]
        p['rectangle_mm']=[float(rect[0]),float(rect[1]),float(rect[2]-rect[0]),float(rect[3]-rect[1])]
        p['letter_slot_mm']=[float(rect[0]),float(rect[1]),LETTER,LETTER]
        # The letter occupies a left corner, not a horizontal row above the content.
        letter=r.fig.text((rect[0]+.1)/160,1-(rect[1]+.15)/r.height,p['panel'][-1],fontsize=11,fontweight='bold',va='top')
        r.strict_text_artists.append((p['panel'],letter))
        # Figure-level source-encoding note belongs to its first panel, inside the boundary.
        notes=[]
        if p['panel'] in ['D01A','D02A']:
            for t in list(r.fig.texts):
                if t.get_text().startswith('Color  relative expression'):
                    t.set_position(((rect[0]+LETTER+PAD)/160,1-(rect[1]+PAD)/r.height));notes.append(t)
        keys,keyheight,names=add_key(r,axes,rect,4.4 if notes else 0);r.panel_keys.extend((p['panel'],lg) for lg in keys)
        r.strict_text_artists.extend((p['panel'],t) for lg in keys for t in lg.get_texts())
        if notes:keyheight=max(keyheight,4.4);r.strict_text_artists.extend((p['panel'],t) for t in notes)
        r.content_headers[p['panel']]=(keys,notes)
        r.legend_ownership.append(dict(panel=p['panel'],dataset_names=names,legend_rows=[dict(kind='dataset' if lg._dataset_key else 'encoding',labels=[getattr(t,'dataset_original',t.get_text()) for t in lg.get_texts()]) for lg in keys]))
        # Component topology is maintained within its authoritative panel rectangle.
        ms=[m for m in r.axes_meta if m['axis'] in axes]
        content=np.array([rect[0]+LETTER,rect[1]+max(keyheight,.2),rect[2]-PAD,rect[3]-PAD])
        boxes={str(i):np.array([m['component_rect'][0],m['component_rect'][1],m['component_rect'][0]+m['component_rect'][2],m['component_rect'][1]+m['component_rect'][3]]) for i,m in enumerate(ms)}
        slots=pack(boxes,content,INTERNAL_GAP) if len(boxes)>1 else {'0':content}
        for i,m in enumerate(ms):
            a=m['axis'];slot=slots[str(i)];m['strict_slot']=slot;m['strict_panel']=p['panel']
            # Out-of-axis identities get a compact, consistent anchor.
            for t in a.texts:
                if t.get_transform()==a.transAxes and t.get_position()[1]>1:
                    xx,_=t.get_position();t.set_position((xx,1.025));t.set_va('bottom')
                    if t.get_text().startswith('• |Δ|'):
                        from matplotlib.transforms import ScaledTranslation
                        t.set_position((xx,1));t.set_transform(a.transAxes+ScaledTranslation(0,26/72,r.fig.dpi_scale_trans))
            if a.axison:
                a.xaxis.labelpad=3;a.yaxis.labelpad=3
                if a.get_xlabel() and not hasattr(a.xaxis.label,'dataset_slots'):a.set_xlabel(r.wrap(a.get_xlabel(),max(12,slot[2]-slot[0]-2)))
                if a.get_ylabel() and not hasattr(a.yaxis.label,'dataset_slots'):a.set_ylabel(r.wrap(a.get_ylabel(),max(13,(slot[3]-slot[1])*.8)))
                # First fit uses generous but explicit text margins; subsequent passes use measured extents.
                p0=axisbox(r,a);bb=a.get_tightbbox(r.fig.canvas.get_renderer()).transformed(r.fig.dpi_scale_trans.inverted());tight=np.array([bb.x0*25.4,r.height-bb.y1*25.4,bb.x1*25.4,r.height-bb.y0*25.4]);margins=np.maximum([p0[0]-tight[0],p0[1]-tight[1],tight[2]-p0[2],tight[3]-p0[3]],0)
                margins=np.minimum(margins,[max(8,(slot[2]-slot[0])*.55),10,12,13])
                pos=slot+np.array([margins[0],margins[1],-margins[2],-margins[3]])
                if pos[2]-pos[0]<4:pos[0]=slot[0]+min(8,(slot[2]-slot[0])*.45);pos[2]=slot[2]-1
                if pos[3]-pos[1]<4:pos[1]=slot[1]+1;pos[3]=slot[3]-5
            else:pos=slot.copy()
            setbox(r,a,pos)
        r.strict_groups.append((p,axes,rect,keys))
    for _ in range(10):
        r.fig.canvas.draw();ren=r.fig.canvas.get_renderer()
        for m in r.axes_meta:
            a=m['axis'];slot=m['strict_slot'];p0=axisbox(r,a)
            b=axis_contentbox(r,a,ren)
            for parent,cb in r.colorbars:
                if parent is a:
                    b[3]=max(b[3],p0[3]+7.1)
            if a.get_aspect()!='auto':
                d=np.maximum([slot[0]-b[0],slot[1]-b[1],b[2]-slot[2],b[3]-slot[3]],0)
            else:d=np.array([slot[0]+PAD-b[0],slot[1]+PAD-b[1],b[2]-slot[2]+PAD,b[3]-slot[3]+PAD])
            pp=p0+np.array([d[0],d[1],-d[2],-d[3]])*.7
            if pp[2]-pp[0]<3 or pp[3]-pp[1]<2:raise RuntimeError((r.fid,m['strict_panel'],'text cannot fit',slot,b,pp))
            setbox(r,a,pp)
        sync_edges(r)
    sync_edges(r)
    # Equal-aspect embeddings use less width after height alignment. Their shared-y
    # marginal distribution starts at the actual embedding edge, not its old slot.
    for link in r.y_links:
        aa=link['axes']
        if len(aa)==2 and aa[0].get_aspect()!='auto' and not aa[0].axison:
            left=axisbox(r,aa[0]);right=axisbox(r,aa[1]);right[0]=left[2]+INTERNAL_GAP;setbox(r,aa[1],right)
    if r.fid in ['D01','D02']:
        # Marker columns and the library strip refer to exactly the same states.
        aa=[m['axis'] for m in r.axes_meta if m['panel']==r.fid+'A'];bb=[axisbox(r,a) for a in aa]
        lo=max(b[0] for b in bb);hi=min(b[2] for b in bb)
        for a,b in zip(aa,bb):b[0]=lo;b[2]=hi;setbox(r,a,b)
    fit_content_envelopes(r)
    # Thin numeric tick labels only when their measured glyph boxes collide.
    r.fig.canvas.draw();ren=r.fig.canvas.get_renderer()
    for m in r.axes_meta:
        a=m['axis']
        if not a.axison:continue
        for which in ['x','y']:
            axis=a.xaxis if which=='x' else a.yaxis
            ticks=a.get_xticks() if which=='x' else a.get_yticks();labels=axis.get_ticklabels()
            try:[float(t.get_text().replace('−','-')) for t in labels if t.get_text()]
            except ValueError:continue
            bb=[t.get_window_extent(ren) for t in labels if t.get_visible() and t.get_text()]
            d=0 if which=='x' else 1
            bb.sort(key=lambda b:b.x0 if d==0 else b.y0)
            if which=='x' and any(t.get_rotation()%180 for t in labels):continue
            if len(bb)>2 and any((bb[i].x1>bb[i+1].x0-2 if d==0 else bb[i].y1>bb[i+1].y0-2) for i in range(len(bb)-1)):
                keep=ticks[::2]
                if len(keep)>4:keep=ticks[::max(2,int(np.ceil(len(ticks)/4)))]
                (a.set_xticks if which=='x' else a.set_yticks)(keep)
    fit_visible_edges(r)
    r.fig.canvas.draw();ren=r.fig.canvas.get_renderer();violations=[];boxes=[]
    record_content_edges(r,ren)
    for p,axes,rect,keys in r.strict_groups:
        artists=[]
        for a in axes:artists+=text_items(a)
        artists += [t for pid,t in r.strict_text_artists if pid==p['panel']]
        for parent,cb in r.colorbars:
            if parent in axes:artists+=text_items(cb)
        for t in dict.fromkeys(artists):
            bb=extent(r,t,ren);text=getattr(t,'dataset_original',t.get_text())
            if not np.isfinite(bb).all():continue
            boxes.append(dict(panel=p['panel'],text=text,fontsize=t.get_fontsize(),bbox_mm=bb.tolist()))
            excess=max(rect[0]-bb[0],rect[1]-bb[1],bb[2]-rect[2],bb[3]-rect[3])
            if excess>.15:violations.append(dict(panel=p['panel'],text=text,excess_mm=float(excess),bbox_mm=bb.tolist()))
        for lg in keys:
            bb=extent(r,lg,ren);excess=max(rect[0]-bb[0],rect[1]-bb[1],bb[2]-rect[2],bb[3]-rect[3])
            if excess>.15:violations.append(dict(panel=p['panel'],legend=True,excess_mm=float(excess)))
    r.strict_textboxes=boxes;r.strict_violations=violations
    return [p for p,aa,rect,keys in r.strict_groups]

def sync_edges(r):
    # Physical plot edges are shared only by components encoding the same y variable.
    for link in r.y_links:
        axes=link['axes'];bb=[axisbox(r,a) for a in axes]
        top=max(b[1] for b in bb);bottom=min(b[3] for b in bb)
        for a,b in zip(axes,bb):b[1]=top;b[3]=bottom;setbox(r,a,b)
    for parent,cb in r.colorbars:
        p=axisbox(r,parent);old=axisbox(r,cb);setbox(r,cb,[p[0],p[3]+2.8,min(p[0]+19,p[2]),p[3]+4.2])

def fit_content_envelopes(r):
    # Fixed-size legends/text stay fixed-size; only axis positions and spans change.
    for _ in range(24):
        r.fig.canvas.draw();ren=r.fig.canvas.get_renderer()
        worst=0
        for p,axes,rect,keys in r.strict_groups:
            keys,notes=r.content_headers[p['panel']]
            header=list(keys)+list(notes)
            if header:
                hb=unionbox([extent(r,t,ren) for t in header]);dy=rect[1]+PAD-hb[1]
                for lg in keys:
                    anchor=lg.get_bbox_to_anchor().transformed(r.fig.transFigure.inverted())
                    lg.set_bbox_to_anchor((anchor.x0,anchor.y0-dy/r.height),transform=r.fig.transFigure)
                for t in notes:
                    x,y=t.get_position();t.set_position((x,y-dy/r.height))
                top=hb[3]+dy+1.2
            else:top=rect[1]+PAD
            bottom=rect[3]-PAD
            allaxes=list(axes)+[cb for parent,cb in r.colorbars if parent in axes]
            box=unionbox([axis_contentbox(r,a,ren) for a in allaxes])
            worst=max(worst,abs(box[1]-top),abs(box[3]-bottom))
            scale=(bottom-top)/(box[3]-box[1])
            for a in axes:
                b=axisbox(r,a);b[[1,3]]=top+(b[[1,3]]-box[1])*scale
                if b[3]-b[1]<2:raise RuntimeError((p['panel'],'full content cannot fit',box,rect))
                setbox(r,a,b)
        sync_edges(r)
        if worst<.008:break

def record_content_edges(r,ren):
    r.content_envelopes=[]
    visible=visible_content_boxes(r)
    for p,axes,rect,keys in r.strict_groups:
        allaxes=list(axes)+[cb for parent,cb in r.colorbars if parent in axes]
        keys,notes=r.content_headers[p['panel']]
        box=unionbox([axis_contentbox(r,a,ren) for a in allaxes]+[extent(r,t,ren) for t in list(keys)+list(notes)])
        r.content_envelopes.append(dict(panel=p['panel'],rectangle_mm=p['rectangle_mm'],content_bbox_mm=visible[p['panel']].tolist(),artist_bbox_mm=box.tolist(),definition='Rendered non-background content at 300 dpi excluding only the panel-letter slot; axes, ticks, labels, legends, annotations and colorbars included'))
    (r.out/'verification'/f'{r.fid}.content-envelopes.json').write_text(json.dumps(r.content_envelopes,ensure_ascii=False,indent=2))
    r.row_align_checks=[]
    for edge in [1,3]:
        groups={}
        for rec in r.content_envelopes:
            x,y,w,h=rec['rectangle_mm'];yy=y if edge==1 else y+h
            groups.setdefault(round(yy,4),[]).append(rec)
        for yy,values in groups.items():
            if len(values)<2:continue
            actual=[v['content_bbox_mm'][edge] for v in values]
            r.row_align_checks.append(dict(edge='top' if edge==1 else 'bottom',panel_edge_mm=yy,panels=[v['panel'] for v in values],measured_edges_mm=actual,max_y_edge_difference_mm=float(max(actual)-min(actual)),basis='Complete content envelope, excluding letters'))

def visible_content_boxes(r):
    pixels=np.asarray(r.fig.canvas.buffer_rgba())[:,:,:3]
    scale=r.fig.dpi/25.4;result={}
    for p,axes,rect,keys in r.strict_groups:
        x0,y0,x1,y1=np.round(rect*scale).astype(int)
        mask=np.min(pixels[y0:y1,x0:x1],axis=2)<245
        # The reserved letter corner contains no scientific content.
        mask[:int(np.ceil(LETTER*scale)),:int(np.ceil(LETTER*scale))]=False
        yy,xx=np.where(mask)
        assert len(xx),(p['panel'],'empty content')
        result[p['panel']]=np.array([x0+xx.min(),y0+yy.min(),x0+xx.max()+1,y0+yy.max()+1])/scale
    return result

def fit_visible_edges(r):
    # Glyph ink does not fill its font bbox. Align the rendered content itself too.
    r.fig.set_dpi(300)
    for _ in range(12):
        r.fig.canvas.draw();ren=r.fig.canvas.get_renderer();visible=visible_content_boxes(r);worst=0
        for p,axes,rect,keys in r.strict_groups:
            b=visible[p['panel']];top=rect[1]+CONTENT_PAD;bottom=rect[3]-CONTENT_PAD
            dt=top-b[1];db=bottom-b[3];worst=max(worst,abs(dt),abs(db))
            keys,notes=r.content_headers[p['panel']]
            if keys or notes:
                for lg in keys:
                    anchor=lg.get_bbox_to_anchor().transformed(r.fig.transFigure.inverted())
                    lg.set_bbox_to_anchor((anchor.x0,anchor.y0-dt/r.height),transform=r.fig.transFigure)
                for t in notes:
                    x,y=t.get_position();t.set_position((x,y-dt/r.height))
                allaxes=list(axes)+[cb for parent,cb in r.colorbars if parent in axes]
                inner=unionbox([axis_contentbox(r,a,ren) for a in allaxes]);origin=inner[1]
                factor=(inner[3]-origin+db)/(inner[3]-origin)
                for a in axes:
                    p0=axisbox(r,a);p0[[1,3]]=origin+(p0[[1,3]]-origin)*factor;setbox(r,a,p0)
            else:
                factor=(bottom-top)/(b[3]-b[1])
                for a in axes:
                    p0=axisbox(r,a);p0[[1,3]]=top+(p0[[1,3]]-b[1])*factor;setbox(r,a,p0)
        sync_edges(r)
        if worst<.085:break
