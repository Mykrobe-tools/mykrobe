class CreateAndSaveMixin(object):

    @classmethod
    def create_and_save(cls, *args, **kwargs):
        return cls.create(*args, **kwargs).save()
