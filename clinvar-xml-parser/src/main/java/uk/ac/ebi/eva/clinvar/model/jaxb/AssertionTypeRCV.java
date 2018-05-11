//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2018.05.09 at 11:38:04 AM BST 
//


package uk.ac.ebi.eva.clinvar.model.jaxb;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;


/**
 * Assertion is used to represent the type of relationship between the trait set and the measure set. This is stored in
 * 						GTR.clinvar.measure_trait.relat_type.
 * 			
 * 
 * <p>Java class for AssertionTypeRCV complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="AssertionTypeRCV">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="Type" use="required" type="{}AssertionTypeAttr" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "AssertionTypeRCV")
public class AssertionTypeRCV {

    @XmlAttribute(name = "Type", required = true)
    protected AssertionTypeAttr type;

    /**
     * Gets the value of the type property.
     * 
     * @return
     *     possible object is
     *     {@link AssertionTypeAttr }
     *     
     */
    public AssertionTypeAttr getType() {
        return type;
    }

    /**
     * Sets the value of the type property.
     * 
     * @param value
     *     allowed object is
     *     {@link AssertionTypeAttr }
     *     
     */
    public void setType(AssertionTypeAttr value) {
        this.type = value;
    }

}
