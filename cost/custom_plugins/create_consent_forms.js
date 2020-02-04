// create html for consent forms, and questions with required consent to proceed
// create consent form, data consent form, load the declaration of consent for data protections
 
function create_consent_form(){
  consent_text = "<p style='font-size:30px'>Informed Consent Form for Experiments on Cognitive Control</p>"
  consent_text += "<p style='font-size:25px'>Kyb Lab at the Max Planck Institute for Biological Cybernetics</p>"
  consent_text += "<p style='font-size:20px'>This is a psychology experiment being conducted by Dr. Peter Dayan, director of the Max Planck Institute for Biological Cybernetics, and the members of his lab. In order to consent to participate, you MUST meet the following criteria: </p>"
  consent_text += "<p style='font-size:20px'>-   18 years of age or older. </p>"
  consent_text += "<p style='font-size:20px'>-   Fluent speaker of English. </p>"
  consent_text += "<p style='font-size:20px'>-   Have not previously participated in this experiment. </p>"
  consent_text += "<p style='font-size:20px'>This study is designed to look at how people learn how to make decisions to accomplish their goals. In this task, you will be asked to make choices, play games, and answer questions related to those games.  The study will take about 30 minutes and will pay $2 plus a performance-dependent bonus of up to a maximum of $10. The performance bonus is explained in more detail in the instructions that follow. </p>"
  consent_text += "<p style='font-size:20px'>Your participation in this research is voluntary. You may refrain from answering any questions that make you uncomfortable and may withdraw your participation at any time without penalty by exiting this task and alerting the experimenter. You may choose not to complete certain parts of the task or answer certain questions. You may contact us at the address provided below if you have additional questions or concerns.</p>"
  consent_text += "<p style='font-size:20px'>Other than monetary compensation, participating in this study will provide no direct benefits to you. But we hope that this research will benefit society at large by contributing towards establishing a scientific foundation for improving human learning and cognitive control abilities.</p>"
  consent_text += "<p style='font-size:20px'>Your online username may be connected to your individual responses, but we will not be asking for any additional personally identifying information, and we will handle responses as confidentially as possible. We cannot however guarantee the confidentiality of information transmitted over the Internet. We will be keeping deidentified data collected as part of this experiment indefinitely. Data used in scientific publications will remain completely anonymous. </p>"
  consent_text += "<p style='font-size:20px'>If you have any questions about the study, feel free to contact our lab. Dr. Dayan and his lab members can be reached at kyblab.tuebingen@gmail.com. </p>"
  consent_text += "<p style='font-size:20px'>By selecting the consent option below, I acknowledge that I am 18 or older, that I am a fluent speaker of English, that I have read this consent form, and that I agree to take part in the research. </p>"
  return consent_text
}

function create_data_consent_form(){
  consent_text = "<p style='font-size:30px'>Legal Declaration of Consent under Data Protection Law</p>"
  consent_text = "<p style='font-size:25px'>Online Experiments on Cognitive Control</p>"
  consent_text += "<a href='Datenschutz_online.pdf' target='_blank'> Please click here to view our Data Protection Information Sheet, yours to keep.</a>"
  consent_text += "<p>I have received and took note of the written <strong>Data Protection Information Sheet</strong> for this study. In doing so, I had sufficient time and opportunity to ask questions about data protection and reconsider my participation in the study.</p>"
  consent_text += "<p></p>"
  consent_text += "<p>I am aware that:</p>"
  consent_text += "<p>- the processing and use of the collected data occurs in a pseudoanonymised form within the scope of the legally prescribed provisions. As a general rule, the storage occurs in the form of answered questionnaires, as well as electronic data, for a duration of 10 years or longer, if this is required by the purpose of the study.</p>"
  consent_text += "<p>- by providing of further personal data in pseudoanonymised form, collected personal data may be used for the preparation of anonymised scientific research work and may also be published and used in an anonymised form in medical journals and scientific publications, so that a direct assignment to my person cannot be established.</p>"
  consent_text += "<p>- the information obtained during the course of this study may also be sent in an anonymised form to cooperation partners within the scope of the European General Data Protection Regulation for scientific purposes and to cooperation partners outside of the European Union, i.e. to countries with a lower data protection level (this also applies to the USA).</p>"
  consent_text += "<p>- the data collected within the scope of the study can also be used and processed in the future inside of the Max Planck Institute.</p>"
  consent_text += "<p></p>"
  consent_text += "<p>I was informed about my rights, that at any time:</p>"
  consent_text += "<p>- I can withdraw this declaration of consent.</p>"
  consent_text += "<p>- I can request information about my stored data and request the correction or blocking of data.</p>"
  consent_text += "<p>- By cancellation of my participation in the study, I can request that any personal data of mine collected until that point are immediately deleted or anonymised.</p>"
  consent_text += "<p>- I can request that my personal data are handed out to me or to third parties (if technically feasible).</p>"
  consent_text += "<p></p>"
  consent_text += "<p>I herewith declare that:</p>"
  consent_text += "<p>- I have been adequately informed about the collection and processing of my personal data and rights.</p>"
  consent_text += "<p>- I consent to the collection and processing of personal data within the scope of the study and its pseudoanonymised disclosure, so that only the persons conducting the study can establish a link between the data and my person.</p>"
  return consent_text
};

function create_consent_timeline(){ //make timelines with non-optional responses
  var consent_form = {
    type: 'survey-multi-select',
    preamble: create_consent_form(),
    required_message: "You must consent to start the experiment.",
    button_label: "Continue",
    questions: [
      {
        prompt: "", 
        options: ["I consent"], 
        horizontal: false,
        required: true,
        name: "consent"
      }
    ], 
  };

  var data_protections = {
    type: 'survey-multi-select',
    preamble: create_data_consent_form(),
    required_message: "You must consent to start the experiment.",
    button_label: "Start Experiment",
    questions: [
    {
      prompt: "",
      options: ["I agree to participate in this experiment"],
      required: true,
      name: "data consent"
    },
    {
      prompt: "",
      options: ["I consent to the use of my data described in the Data Protection Information Sheet and confirm having received a copy of the Data Protection Sheet"],
      required: true,
      name: "data info sheet"
    },
    {
      prompt: "",
      options: ["I consent to data transfer from the MPI for Biological Cybernetics encrypted database to the project-related collaborators: inside of the Max Planck Society and affiliated research institutes, or at partnering institutions like the University of Tuebingen and New York University"],
      required: true,
      name: "collaborators"
    }]
  };

  return timeline = [consent_form, data_protections]

}