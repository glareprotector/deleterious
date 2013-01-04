from param import param

import pdb, os, subprocess
import cPickle as pickle
from helper import print_stuff_dec, write_mat, write_vect, super_shorten
import helper
from Bio import AlignIO

import caches
from wrapper_decorator import *
#from wc import wc

# since wrappers and caches are objects that can be cached, should record the wrapper that is constructing them.

# note for all caches: recalculate specifies whether you want to use the files/pickles that are present when cache is created

# act of caching: moving file from temporary holding spot to proper folder.  cache structure: file system

# this is the function wrapper

class wrapper(object):

    def get_max_cache_size(self):
        return 5

    def get_to_pickle(self):
        return True

    @classmethod
    def get_all_keys(cls, params, self = None):
        return set()

    @classmethod
    def get_object_key(cls, params, self = None):
        import param
        import pdb

        return params.restriction(params, cls.get_all_keys(params, self))

    def whether_to_override(self, object_key):
        return True

    def super_shorten(self):
        return True

    def specificity(self):
        return None

    def is_indexed(self):
        return True

    def makes_index(self):
        return False

    def __repr__(self):
        return self.__class__.__name__

    def get_folder(self, object_key):
        return global_stuff.BIN_FOLDER

    # if the wrapper index stuff, its objects can obtain indicies thru maker->object_key_to_index
    # in general, the methods i make should have a maker pointer so that params can be accessed thru the maker
    def get_name(self, object_key, to_reindex = global_stuff.to_reindex):
        object_key = self.get_object_key(object_key)
        if not self.super_shorten():
            assert False
            if not self.makes_index():
                return helper.shorten(str(object_key))
            elif not to_reindex:
                return helper.shorten(str(object_key))
            else:
                #pdb.set_trace()
                assert self.maker.object_key_to_index.has(object_key, to_reindex) == True
                return str(self.maker.object_key_to_index.get(object_key))
        else:
            if not self.makes_index():
                return helper.super_shorten(str(object_key))
            elif not to_reindex:
                return helper.super_shorten(str(object_key))
            else:
                #pdb.set_trace()
                assert self.maker.object_key_to_index.has(object_key, to_reindex) == True
                return global_stuff.super_shorten(str(self.maker.object_key_to_index.get(object_key)))
            

    # can't get and set same parameter in the same run of constructor.  at the beginning of constructor, set maker.  in init of whatever object, would fetch maker.  can NOT get, then set.  however, can set and THEN get.
    def basic_init(self, maker, params):
        self.temp_used_keys = []
        self.temp_dependents_keys = []
        self.temp_new_param_keys = []
        self.temp_dependents_all_keys_key_key_sets = []
        maker.set_param(params, "source_instance", self)
        self.used_keys_cache = caches.ukcO(maker, params)
        self.set_keys_cache = caches.skcO(maker, params)
        self.all_keys_key_key_set_cache = caches.akkkscO(maker, params)
        
        self.all_keys_cache = caches.akcO(maker, params)
        if self.makes_index():
            self.object_key_to_index = caches.index_cache(maker, params)
        self.maker = maker

    def before_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        return maker, params
    
# if set param, means its value was based on other params, so only other params matter.
    def __init__(self, maker, params = param({})):
#        pdb.set_trace()
        maker, params = self.before_init(maker, params)
        self.base_folder = global_stuff.base_folder
        self.other_init(maker, params)
        self.basic_init(maker, params)

    def other_init(self, maker, params):
        pass

    def set_param(self, params, key, val):
        params.set_param(key, val)
        self.temp_new_param_keys[-1].add(key)
        return params

    def get_param(self, params, key, record = True):
        if record:
            self.temp_used_keys[-1].add(key)
        return params.get_param(key)

    # wrapper refers to which class, not a specific wrapper instance
    @print_stuff_dec
    def get_var_or_file(self, wrapper, params, recalculate = False, to_pickle = False, to_filelize = False, always_recalculate = False):
#        print '                      called  ', self, wrapper, recalculate, to_pickle, to_filelize, always_recalculate
        #pdb.set_trace()
        self.set_param(params, "which_wrapper_class", wrapper)
        from wc import wc
#        print wc.object_key_to_index.dump
        the_wrapper = self.old_get_var_or_file(wc, params, True, False, False)
        used_keys, all_keys, x, all_keys_key_key_set = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize, always_recalculate)
        self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        self.temp_dependents_all_keys_key_key_sets[-1] = self.temp_dependents_all_keys_key_key_sets[-1].union(all_keys_key_key_set)
#        print '                      returned ', self, wrapper, recalculate, to_pickle, to_filelize, always_recalculate
        return x

    def old_get_var_or_file(self, the_wrapper, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False):

        used_keys, all_keys, x, all_keys_key_key_set = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize, always_recalculate)

        self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        self.temp_dependents_all_keys_key_key_sets[-1] = self.temp_dependents_all_keys_key_key_sets[-1].union(all_keys_key_key_set)
        return x




    # returns only the object, but after decorating, will return used_keys, all_keys, object
    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        pass

    def get_cache(self):
        pass

    def has(self, object_key, recalculate, check_remote = False):
        return self.cache.has(object_key, recalculate, check_remote)

    def get(self, object_key, recalculate):
        return self.cache.get(object_key, recalculate)

    def set(self, object_key, object, to_pickle, params):
        return self.cache.set(object_key, object, to_pickle, params)

    def get_file_name(self, object_key):
        location = self.get_file_location(object_key)
        return location.split('/')[-1]

    def record_transferred_file(self, object_key):
        try:
            if not self.was_transferred(object_key):

                f = open(self.get_transferred_files_location(object_key), 'a')
                f.write(self.get_file_name(object_key) + '\n')
                f.close()
        except:
            f = open(self.get_transferred_files_location(object_key), 'w')
            f.write(self.get_file_name(object_key) + '\n')
            f.close()
        

    def was_transferred(self, object_key):
        try:

            folder_transferred_files = helper.get_file_string_set(self.get_transferred_files_location(object_key))
            ans = self.get_file_name(object_key) in folder_transferred_files
            if ans:
                print 'checked remote and found: ', self, object_key
            return ans    
        except:

            return False

    def get_transferred_files_location(self, object_key):
        return self.get_folder(object_key) + 'moved_files'

class shorten_name_wrapper(wrapper):

    def super_shorten(self):
        return True


class always_recalculate_wrapper(wrapper):

    def always_recalculate(self):
        return True
    

class by_pdb_folder_wrapper(wrapper):

    def get_folder(self, object_key):
        return global_stuff.BIN_FOLDER  + object_key.get_param('p') + '_' + object_key.get_param('c') + '_' + str(object_key.get_param('st')) + '_' + str(object_key.get_param('en')) + '/'

    def specificity(self):
        return 'chain'

class experiment_results_wrapper(wrapper):

    def get_folder(self, object_key):
        #return constants.BIN_FOLDER + 'experiment_results/'
        return global_stuff.RESULTS_FOLDER



class indexing_wrapper(wrapper):

    def makes_index(self):
        return True

class indexed_wrapper(wrapper):

    def process_index(self, index):
        self.name = str(index)

class obj_wrapper(wrapper):

    def get_pickle_dumper(self, maker, params):
        return pkdW(maker, params)

    def whether_to_override(self, object_key):
        return True

    def get_file_dumper(self, maker, params):
        return dfdW(maker, params)

    def get_file_location(self, object_key):
        return self.cache.pickle_dumper_wrapper.get_file_location(object_key)
        return self.get_folder(object_key) + self.get_name(object_key) + '.pk'

    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.object_cache_for_wrapper(maker, params)

class mat_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return mfdW(maker, params)
    
    def get_file_dumper(self, maker, params):

        return mfdW(maker, params)
        maker.set_param(params, "which_wrapper_class", wrapper.mfdW)
        return maker.get_var_or_file(wrapper_catalog, params, True, False, False)

class vect_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return vfdW(maker, params)

    def get_file_dumper(self, maker, params):
        return vfdW(maker, params)

class edge_to_int_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return etiW(maker, params)

    def get_file_dumper(self, maker, params):
        return etiW(maker, params)

class int_mat_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return imfdW(maker, params)

    def get_file_dumper(self, maker, params):
        return imfdW(maker, params)

class int_float_tuple_mat_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return iftmfdW(maker, params)

    def get_file_dumper(self, maker, params):
        return iftmfdW(maker, params)

class msa_obj_wrapper(obj_wrapper):

    def get_pickle_dumper(self, maker, params):
        return dadW(maker, params)

    def get_file_dumper(self, maker, params):
        return dadW(maker, params)

class file_wrapper(wrapper):

    def get_backup_location(self):
        return global_stuff.BIN_FOLDER + str(id(self)) + '.backup'

    def get_holding_location(self):
        return global_stuff.get_holding_folder() + str(id(self))

    def get_file_location(self, object_key):

        return self.get_folder(object_key) + self.__repr__() + '=' + self.get_name(object_key)

    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.file_cache_for_wrapper(maker, params)
        
    @print_stuff_dec
    def get_holding_folder(self):
        return global_stuff.get_holding_folder()

    
# for now, assume source_wrapper instance is available.  if instance is available, that means it wasn't created yet, so wouldn't be in wrapper registry, so don't need to fetch it from there.  in the future, might separately create source instance first, in which case might instead pass in the parameters needed to fetch the param.  or could fetch the instance first externally, then pass to this __init__
class generic_dumper_wrapper(file_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):

        return self.source_wrapper.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return self.source_wrapper.whether_to_override(object_key)

    def __repr__(self):
        return self.__class__.__name__  + '*' + self.source_wrapper.__repr__()


    def get_name(self, object_key, to_reindex = global_stuff.to_reindex):

        return self.source_wrapper.get_name(object_key)
        
    def other_init(self, maker, params):
        self.source_wrapper = maker.get_param(params, "dumper_source_instance")
        #print self, self.source_wrapper
        #pdb.set_trace()
        maker.set_param(params, "source_instance", self)
        self.cache = caches.file_cache_for_wrapper(maker, params)

    def before_init(self, maker, params):
        return maker, params

    def dump_object(self, obj):
        pass

    def get_folder(self, object_key):
        return self.source_wrapper.get_folder(object_key)

    # never pickle, since object returned is a file handle
    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        assert to_pickle == False
        #assert recalculate == True
        # manually set always_recalculate to false here
        # always_recalculate = False
        # if always recalculate is true, have option here to set always_recalculate to be true.  where to decide this?  could be property of source_wrapper.
        # if you call dumper directly, if always_recalculate was already true, then recalculate should be true here.  if was part of caching, and always_recalculate is True, it should be changed to false here
        if always_recalculate:
            if always_recalculate != 2:
                always_recalculate = False
        #ERROR?
        obj = self.old_get_var_or_file(self.source_wrapper, params, recalculate, to_pickle, to_filelize, always_recalculate)
        import datetime
        start = datetime.datetime.now()
        self.dump_object(obj)
        #print 'pickle: ', datetime.datetime.now() - start, self.source_wrapper
        return open(self.get_holding_location(), 'rb')

class pkdW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return pickle.load(f)
    
    def dump_object(self, obj):
        pickle.dump(obj, open(self.get_holding_location(), 'wb'))

# assume mat_obj_wrapper reads stuff to float
class mfdW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return helper.read_mat_to_float(f)

    # do i have to create the folder first?
    def dump_object(self, obj):
        write_mat(obj, self.get_holding_location())

class etiW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return helper.read_edge_to_int(f)

    def dump_object(self, obj):
        return helper.write_edge_to_int(obj, self.get_holding_location())

class imfdW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return helper.read_mat_to_int(f)

    def dump_object(self, obj):
        helper.write_mat_raw(obj, self.get_holding_location())

class iftmfdW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return helper.read_mat_to_int_float_tuple(f)

    def dump_object(self, obj):
        helper.write_int_float_tuple_mat(obj, self.get_holding_location())

# assumes vect_obj_wrapper reads stuff to float
class vfdW(generic_dumper_wrapper):

    def read_object_from_file(self, f):
        return helper.read_vect_to_float(f)

    # do i have to create the folder first?
    def dump_object(self, obj):
        write_vect(obj, self.get_holding_location())

class dfdW(generic_dumper_wrapper):

    def read_project_from_file(self, f):
        assert(False)

    def dump_object(self, obj):
        f = open(self.get_holding_location(), 'w')
        f.write(str(obj))
        #pdb.set_trace()
        
class dadW(generic_dumper_wrapper):

    def dump_object(self, object):
        AlignIO.write(object, open(self.get_holding_location(),'w'), 'fasta')

class wrapper_catalog(obj_wrapper):

    def get_max_cache_size(self):
        return 100

    def get_to_pickle(self):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['which_wrapper_class'])

    # since there is no maker, hackishly set it to self
    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.object_cache_for_wrapper(maker, params)

    def is_indexed(self):
        return False

    # params contains which_wrapper, and if which_wrapper is a generic_dumper_wrapper, contains 
    @dec
    def constructor(self, params, recalculate = False, to_pickle = False, to_filelize = False, always_recalculate = False):
#        pdb.set_trace()
        wrapper_instance = self.get_param(params, "which_wrapper_class")(self, params)
        return wrapper_instance

# only purpose of this class is to give wrapper_catalog instance a "maker"
class famished_wrapper(object):

    

    def set_param(self, params, key, val):
        params.set_param(key, val)
        return params

    def get_param(self, params, key, record = True):
        return params.get_param(key)

    def get_var_or_file(self, wrapper, params, recalculate, to_pickle, to_filelize = False):


        self.set_param(params, "which_wrapper_class", wrapper)
        from wc import wc
        wrapper_used_keys, wrapper_all_keys, the_wrapper = self.old_get_var_or_file(wc, params, True, False, False)
        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        #self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        return x

    def old_get_var_or_file(self, the_wrapper, params, recalculate, to_pickle, to_filelize = False):
        #print '                       ', self, wrapper
        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        return x

class experiment_type_wrapper(wrapper):

    def get_folder(self, object_key):
        folder = str(self.get_param(object_key, 'wif', False)) + '_' + str(self.get_param(object_key, 'wob', False)) + '/'
        return global_stuff.RESULTS_FOLDER + folder + self.__repr__() + '/'

class by_uniprot_id_wrapper(wrapper):
    def get_folder(self, object_key):
        folder = str(self.get_param(object_key, 'uniprot_id', False))
        return global_stuff.base_folder + folder + '/'
