"""
TODO: write docstring.
"""

## std
import os
import re
import sys
import time

## ROOT
import ROOT

_files_opened = dict()
verbosity = 1
i_fig = 1
plotsdir = None

#______________________________________________________________________________
def get(filename, path='', rename=''):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global _files_opened
    if _files_opened.has_key(filename):
        f = _files_opened[filename]
    else:
        f = ROOT.gROOT.GetListOfFiles().FindObject(filename)
        if f:
            print 'ipyhep: WARNING: Weird, %s in GetListOfFiles() but not _files_opened.' % filename
            _files_opened[filename] = f
        else:
            f = ROOT.TFile.Open(filename, 'OPEN')
            _files_opened[filename] = f
    if f:
        if path:
            obj = f.Get(path)
            if isinstance(obj, ROOT.TH1):
                if rename:
                    newobj = obj.Clone(rename)
                else:
                    newobj = obj.Clone('%s_copy' % obj.GetName())
                newobj.SetDirectory(0)
                if verbosity >= 1:
                    print 'ipyhep: reading %s:%s/%s' % (filename, path, newobj.GetName())
                    print newobj
                    sys.stdout.flush()
                return newobj
            elif obj:
                ## not histograms, e.g. TTrees, don't clone
                return newobj
            else:
                print 'ipyhep: WARNING get failed on path: %s' % path
        else:
            return f
    else:
        print 'ipyhep: ERROR get failed to open TFile %s' % filename
    return None

#______________________________________________________________________________
def write(obj, filename, dir=''):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global _files_opened
    f = None
    if _files_opened.has_key(filename):
        f = _files_opened[filename]
    else:
        f = ROOT.gROOT.GetListOfFiles().FindObject(filename)
    if not f:
        f = ROOT.TFile.Open(filename, 'RECREATE')
        _files_opened[filename] = f
    d = f.GetDirectory(dir)
    if not d:
        d = make_root_dir(f, dir)
    d.cd()
    if verbosity >= 1:
        print 'ipyhep: writing %s:%s/%s' % (filename, dir, obj.GetName())
        sys.stdout.flush()
    obj.Write()
    return f

#______________________________________________________________________________
def close_all_files():
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global _files_opened
    for fn, f in _files_opened.iteritems():
        if verbosity >= 1:
            print 'ipyhep: closing %s' % fn
            sys.stdout.flush()
        f.Close()

#______________________________________________________________________________
def make_root_dir(f, dir):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dir.rstrip('/')
    dir_split = dir.split('/')
    lead_dir = dir_split[0]
    sub_dirs = dir_split[1:]

    d = f.GetDirectory(lead_dir)
    if not d:
        d = f.mkdir(lead_dir)
    
    if sub_dirs:
        return make_root_dir(d, '/'.join(sub_dirs))
    else:
        return d

#______________________________________________________________________________
def strip_root_ext(filename, exts=None):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    filename = os.path.basename(filename)
    if exts is None:
        exts = [ 
                '.canv.root',
                '.hist.root',
                '.skim.root',
                '.root',
                ]   
    for ext in exts:
        if filename.endswith(ext):
            return filename[:-1*len(ext)]
    return filename

#______________________________________________________________________________
def walk(top, topdown=True):
    """
    os.path.walk like function for TDirectories.
    Return 4-tuple: (dirpath, dirnames, objnames, top)
        dirpath = '/some/path' for some file_name.root:/some/path
        dirnames = ['list', 'of' 'TDirectory', 'keys']
        objnames = ['list', 'of' 'object', 'keys']
        top = this level's TDirectory
    """
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    assert isinstance(top, ROOT.TDirectory)
    fullpath = top.GetPath()  ## /a/b/c.root:/some/path
#    assert fullpath.startswith('/'), fullpath
    assert fullpath.count(':') == 1, fullpath
    dirpath = fullpath.split(':')[1]
    dirnames = []
    objnames = []
    ## filter names for directories
    for k in top.GetListOfKeys():
        n = k.GetName()
        if k.IsFolder():
            dirnames.append(n)
        else:
            objnames.append(n)
    ## sort
    dirnames.sort()
    objnames.sort()
    ## yield
    if topdown:
        yield dirpath, dirnames, objnames, top
    for dn in dirnames:
        d = top.Get(dn)
        for x in walk(d, topdown):
            yield x
    if not topdown:
        yield dirpath, dirnames, objnames, top


#______________________________________________________________________________
def get_obj_list(tfile):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    obj_list = []
    for dirpath, dirnames, objnames, top in walk(tfile):
        for name in objnames:
            obj_path = os.path.join(dirpath, name)
            obj_list.append(obj_path)
    return obj_list


#______________________________________________________________________________
def remove_strange_characters(fn):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    len_fn = len(fn)
    i = 0
    while i < len_fn:
        ch = fn[i]
        if not re.match(r'[a-zA-Z0-9_.\-]', ch): # \w = [a-zA-Z0-9_]
            fn = fn.replace(ch, '$')
            len_fn = len(fn)
        i += 1
    fn = fn.replace('$', '')
    return fn


#______________________________________________________________________________
def save_figures(canvas, name, ext):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global i_fig
    global plotsdir
    saves = None
    if isinstance(ext, str):
        saves = [ext]
    elif isinstance(ext, list):
        saves = ext 
    words = name.split(':')
    clean_words = list()
    for word in words:
        clean_word = word.strip().strip('.')
        if clean_word.endswith('/1000'):
            clean_word = clean_word[:-5]
        clean_words.append(clean_word)
    clean_name = '_'.join(clean_words)
    clean_name = remove_strange_characters(clean_name)
    if not plotsdir:
        timestamp = time.strftime('%Y_%m_%d_%Hh%M')
        plotsdir = 'plots_%s' % timestamp
        if os.path.exists(plotsdir):
            print 'Removing existing plots directory: %s' % plotsdir
            os.system('rm -rf %s' % plotsdir)
        os.makedirs(plotsdir)
    for s in saves:
        fn = '%s/%02i_%s.%s' % (plotsdir, i_fig, clean_name, s)
        assert not os.path.exists(fn)
        canvas.SaveAs(fn)
    i_fig += 1


