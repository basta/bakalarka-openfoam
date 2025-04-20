class mySocketFunctionObject : public functionObjects::fvMeshFunctionObject
{
    // Private data
    word fieldName_;
    label port_;
    word hostName_;

    // Socket communication handler
    autoPtr<socketCommunication> socketComm_;

    // Fields to read
    volVectorField* fieldPtr_;

public:
    // Runtime type information
    TypeName("mySocketFunctionObject");

    // Constructors
    mySocketFunctionObject
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Member Functions
    virtual bool execute();
    virtual bool write();
    virtual bool read(const dictionary& dict);

    // Create the field if it doesn't exist
    void createFieldIfNeeded();

    // Read field data from socket
    void readFieldFromSocket();
};
